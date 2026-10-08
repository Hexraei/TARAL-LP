// Bounded-variable dual revised simplex. Dual steepest-edge pricing of the leaving row, a
// bound-flipping ratio test with Harris tolerances, dual phase 1 by the box-bounded auxiliary
// problem (free -> [-1,1], lower-only -> [0,1], upper-only -> [-1,0], boxed -> [0,0]), cost
// perturbation against dual degeneracy, cost shifting for tiny dual infeasibilities, product-form
// (eta) basis updates and periodic refactorization with the same SparseLU as the primal.
// Variables follow the primal's convention: 0..n-1 structural, n+i the logical of row i (-e_i).
//
// The dual never claims optimality itself: its basis is handed to the primal simplex, which
// restores the true costs, removes any remaining primal or dual infeasibility and certifies the
// optimum with its own tolerances and feasibility check. Infeasibility is reported only from an
// explicit check of the dual ray against the true bounds.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>

#include "taral.hpp"
#include "timing.hpp"

namespace {
constexpr double kPrimalTol = 1e-9;   // leaving-row infeasibility
constexpr double kDualTol = 1e-9;     // reduced-cost feasibility (Harris tolerance)
constexpr double kPivotTol = 1e-9;    // smallest |alpha_rj| considered in the ratio test
constexpr double kPerturb = 5e-7;     // relative cost perturbation
constexpr size_t kRefactorEvery = 120;
constexpr long kMaxIterations = 50'000'000;
constexpr long kStallIterations = 10000;  // without dual objective progress: hand over to the primal

enum Where : char { kBasic, kLower, kUpper, kZero };

struct Eta {
    int r;
    double pivot;
    std::vector<Entry> col;
};

enum class Outcome { Optimal, Infeasible, TimeLimit, Fallback };

class Dual {
public:
    Dual(const Model& md, const std::vector<double>& clo, const std::vector<double>& cup)
        : md_(md), clo_(clo), cup_(cup), m_(int(md.row_names.size())), n_(int(md.col_names.size())) {}
    Outcome run(double time_limit_s, const std::vector<char>* warm);
    long iterations = 0;
    std::vector<char> basis;
    std::string message;
    double factor_s = 0, first_factor_s = 0;  // seconds in LU factorization (all, and the initial one)
    long factors = 0;
    std::string start_diag;  // warm starts only: what the starting basis looked like
    long phase1_iters = 0;
    std::string timeout_info;
    std::string note;  // warm-start fallback, kept apart from the later status messages

private:
    const Model& md_;
    const std::vector<double>&clo_, &cup_;
    int m_, n_;
    std::vector<double> tlo_, tup_;  // true bounds
    std::vector<double> lo_, up_, c_, x_, d_;
    std::vector<int> head_;
    std::vector<Where> where_;
    std::vector<std::vector<Entry>> rows_;  // row-wise structural part of A: (column, value)
    SparseLU lu_;
    std::vector<Eta> etas_;
    SparseLU::Clock::time_point deadline_ = SparseLU::Clock::time_point::max();
    bool timed_out_ = false;  // a factorization hit deadline_

    template <class F>
    void for_col(int j, F f) const {
        if (j < n_)
            for (const Entry& e : md_.cols[j]) f(e.index, e.value);
        else
            f(j - n_, -1.0);
    }
    double dot_col(int j, const std::vector<double>& y) const {
        double s = 0;
        for_col(j, [&](int i, double v) { s += v * y[i]; });
        return s;
    }
    void ftran(std::vector<double>& rhs, std::vector<double>& out) const {
        lu_.ftran(rhs, out);
        for (const Eta& e : etas_) {
            double xr = out[e.r] / e.pivot;
            if (xr != 0)
                for (const Entry& t : e.col) out[t.index] -= t.value * xr;
            out[e.r] = xr;
        }
    }
    void btran(std::vector<double>& v, std::vector<double>& out) const {
        for (auto it = etas_.rbegin(); it != etas_.rend(); ++it) {
            double s = v[it->r];
            for (const Entry& t : it->col) s -= t.value * v[t.index];
            v[it->r] = s / it->pivot;
        }
        lu_.btran(v, out);
    }
    // Nonbasic position that is dual feasible for reduced cost d under the current bounds.
    void place(int j) {
        bool fl = std::isfinite(lo_[j]), fu = std::isfinite(up_[j]);
        if (fl && fu) where_[j] = (lo_[j] == up_[j] || d_[j] >= 0) ? kLower : kUpper;
        else if (fl) where_[j] = kLower;
        else if (fu) where_[j] = kUpper;
        else where_[j] = kZero;
        x_[j] = where_[j] == kLower ? lo_[j] : where_[j] == kUpper ? up_[j] : 0.0;
    }
    bool factor();
    void compute_primal();
    void compute_dual();
    int repair_dual(bool shift);
    bool proves_infeasible(int r, const std::vector<double>& alpha_row, const std::vector<int>& touched) const;
};

bool Dual::factor() {
    const auto t0 = SparseLU::Clock::now();
    struct Tally {
        Dual& d;
        SparseLU::Clock::time_point t0;
        ~Tally() {
            double s = std::chrono::duration<double>(SparseLU::Clock::now() - t0).count();
            if (d.factors++ == 0) d.first_factor_s = s;
            d.factor_s += s;
        }
    } tally{*this, t0};
    for (int attempt = 0; attempt < 5; ++attempt) {
        std::vector<std::vector<Entry>> cols(m_);
        for (int p = 0; p < m_; ++p) for_col(head_[p], [&](int i, double v) { cols[p].push_back({i, v}); });
        etas_.clear();
        if (lu_.factor(m_, cols, deadline_)) return true;
        if (lu_.timed_out) {
            timed_out_ = true;
            timeout_info = "after " + std::to_string(lu_.pivots_done) + " of " + std::to_string(m_) + " pivots, " +
                           std::to_string(lu_.active_nnz) + " entries left in the active submatrix";
            return false;
        }
        for (int p : lu_.bad_pos) where_[head_[p]] = kLower;  // placed properly after compute_dual
        for (size_t t = 0; t < lu_.bad_pos.size(); ++t) {
            int j = n_ + lu_.bad_rows[t];
            head_[lu_.bad_pos[t]] = j;
            where_[j] = kBasic;
        }
    }
    return false;
}

void Dual::compute_primal() {
    std::vector<double> rhs(m_, 0.0), xb;
    for (int j = 0; j < n_ + m_; ++j)
        if (where_[j] != kBasic && x_[j] != 0) for_col(j, [&](int i, double v) { rhs[i] -= v * x_[j]; });
    ftran(rhs, xb);
    for (int p = 0; p < m_; ++p) x_[head_[p]] = xb[p];
}

void Dual::compute_dual() {
    std::vector<double> cb(m_), y;
    for (int p = 0; p < m_; ++p) cb[p] = c_[head_[p]];
    btran(cb, y);
    for (int j = 0; j < n_ + m_; ++j) d_[j] = where_[j] == kBasic ? 0.0 : c_[j] - dot_col(j, y);
}

// Restore dual feasibility of the nonbasics: boxed variables move to the bound their reduced cost
// asks for; others either have their cost shifted (shift) or are counted as dual infeasible.
// Returns the number of remaining dual infeasibilities. x_B must be recomputed afterwards.
int Dual::repair_dual(bool shift) {
    int bad = 0;
    for (int j = 0; j < n_ + m_; ++j) {
        if (where_[j] == kBasic) continue;
        bool fl = std::isfinite(lo_[j]), fu = std::isfinite(up_[j]);
        if (fl && fu) {
            // Keep a boxed nonbasic at its current bound unless its reduced cost is wrong-signed by more than the
            // dual tolerance. A sign flip from rounding noise would otherwise move it to the other bound (a full-range
            // move), so the next pivot undoes it and the loop refactorizes every iteration.
            if (lo_[j] != up_[j] && ((where_[j] == kLower && d_[j] >= -kDualTol) || (where_[j] == kUpper && d_[j] <= kDualTol))) {
                x_[j] = where_[j] == kLower ? lo_[j] : up_[j];
                continue;
            }
            place(j);
            continue;
        }
        double v = where_[j] == kLower ? std::min(0.0, d_[j]) : where_[j] == kUpper ? std::max(0.0, d_[j]) : d_[j];
        if (!fl && !fu) where_[j] = kZero, x_[j] = 0;
        else place(j);
        if (std::abs(v) <= kDualTol) continue;
        if (shift) c_[j] -= v, d_[j] -= v;
        else ++bad;
    }
    return bad;
}

// Row r of B^{-1}A gives x_r = -sum_{j nonbasic} alpha_rj x_j. If the range of that sum over the
// TRUE bounds misses [lo_r, up_r] by a clear margin, the model is infeasible.
bool Dual::proves_infeasible(int r, const std::vector<double>& alpha_row, const std::vector<int>& touched) const {
    double mn = 0, mx = 0, scale = 0;
    for (int j : touched) {
        if (where_[j] == kBasic) continue;
        double a = -alpha_row[j];
        if (a == 0) continue;
        double l = tlo_[j], u = tup_[j];
        double vmin = a > 0 ? a * l : a * u, vmax = a > 0 ? a * u : a * l;
        mn += vmin, mx += vmax;
        if (std::isfinite(vmin)) scale = std::max(scale, std::abs(vmin));
        if (std::isfinite(vmax)) scale = std::max(scale, std::abs(vmax));
    }
    int jr = head_[r];
    double margin = 1e-6 * (1 + scale);
    return (std::isfinite(tlo_[jr]) && mx < tlo_[jr] - margin) || (std::isfinite(tup_[jr]) && mn > tup_[jr] + margin);
}

Outcome Dual::run(double time_limit_s, const std::vector<char>* warm) {
    Stopwatch elapsed;
    int N = n_ + m_;
    tlo_.assign(N, 0), tup_.assign(N, 0), c_.assign(N, 0), x_.assign(N, 0), d_.assign(N, 0);
    where_.assign(N, kLower), head_.resize(m_);
    for (int j = 0; j < n_; ++j) tlo_[j] = clo_[j], tup_[j] = cup_[j], c_[j] = md_.cost[j];
    for (int i = 0; i < m_; ++i) tlo_[n_ + i] = md_.row_lo[i], tup_[n_ + i] = md_.row_up[i];
    auto keep_basis = [&] { basis.assign(where_.begin(), where_.end()); };
    for (int j = 0; j < N; ++j)
        if (tlo_[j] > tup_[j]) {
            message = "inconsistent bounds";
            return Outcome::Fallback;
        }
    rows_.assign(m_, {});
    for (int j = 0; j < n_; ++j)
        for (const Entry& e : md_.cols[j]) rows_[e.index].push_back({j, e.value});

    lo_ = tlo_, up_ = tup_;
    auto slack_basis = [&] {
        std::fill(where_.begin(), where_.end(), kLower);
        for (int i = 0; i < m_; ++i) head_[i] = n_ + i, where_[n_ + i] = kBasic;
    };
    const auto start = SparseLU::Clock::now();
    const auto limit = std::chrono::duration_cast<SparseLU::Clock::duration>(std::chrono::duration<double>(time_limit_s));
    deadline_ = start + limit;
    bool warmed = warm && int(warm->size()) == N && std::count(warm->begin(), warm->end(), char(kBasic)) == m_;
    if (warmed) {
        for (int j = 0, p = 0; j < N; ++j) {
            where_[j] = Where((*warm)[j]);
            if (where_[j] == kBasic) head_[p++] = j;
        }
        deadline_ = start + limit / 4;  // a warm basis gets a quarter of the limit to factor
    } else {
        slack_basis();
    }
    bool factored = factor();
    if (!factored && warmed && timed_out_) {
        note = "warm basis not factored within " + std::to_string(time_limit_s / 4) + " s (" + timeout_info + "), slack basis used";
        std::fprintf(stderr, "dual: %s\n", note.c_str());
        slack_basis();
        deadline_ = start + limit;
        timed_out_ = false;
        factored = factor();
    } else deadline_ = start + limit;
    if (!factored) {
        if (timed_out_) return Outcome::TimeLimit;
        message = "initial basis could not be factored";
        return Outcome::Fallback;
    }
    // Cost perturbation in the direction that widens the current dual feasibility.
    {
        compute_dual();
        std::mt19937_64 rng(2718281);
        double base = kPerturb;  // also for zero costs, which are otherwise entirely dual degenerate
        for (int j = 0; j < N; ++j) {
            bool fl = std::isfinite(lo_[j]), fu = std::isfinite(up_[j]);
            if (lo_[j] == up_[j] || (!fl && !fu)) continue;
            double xi = base * (1 + std::abs(c_[j])) * (1 + double(rng() % 1000000) / 1e6);
            bool lower = fl && fu ? (where_[j] == kBasic ? c_[j] >= 0 : d_[j] >= 0) : fl;
            c_[j] += lower ? xi : -xi;
        }
    }
    compute_dual();
    bool phase1 = false;
    const int dual_inf0 = repair_dual(false);
    if (dual_inf0 > 0) {  // dual phase 1 on the box-bounded auxiliary problem
        phase1 = true;
        for (int j = 0; j < N; ++j) {
            bool fl = std::isfinite(tlo_[j]), fu = std::isfinite(tup_[j]);
            lo_[j] = fl && fu ? 0.0 : fl ? 0.0 : fu ? -1.0 : -1.0;
            up_[j] = fl && fu ? 0.0 : fl ? 1.0 : fu ? 0.0 : 1.0;
        }
        for (int j = 0; j < N; ++j)
            if (where_[j] != kBasic) place(j);
    }
    compute_primal();
    if (warmed) {
        int pinf = 0;
        double psum = 0;
        for (int p = 0; p < m_; ++p) {
            int j = head_[p];
            double v = std::max({0.0, lo_[j] - x_[j], x_[j] - up_[j]});
            if (v > kPrimalTol) ++pinf, psum += v;
        }
        start_diag = "start: " + std::to_string(dual_inf0) + " dual infeasibilities" + (phase1 ? " (phase 1)" : "") + ", " +
                     std::to_string(pinf) + " primal infeasible basics (sum " + std::to_string(psum) + ")";
    }

    std::vector<double> weight(m_, 1.0);  // dual steepest-edge weights ||e_r' B^{-1}||^2 (by position)
    std::vector<double> alpha_row(N, 0.0), unit, rho, col(m_), alpha, tau, delta(m_);
    std::vector<char> mark(N, 0);
    std::vector<int> touched;
    struct Cand {
        int j;
        double a, ratio, harris;
    };
    std::vector<Cand> cands;
    std::vector<int> flips;
    bool fresh = true;
    int shifts = 0;
    double best_obj = -kInf;  // dual objective (c'x of the current basic solution) at refactorizations
    long progress_at = 0;
    long refreshes = 0;  // refactorizations so far (a storm of them means the row/column computations disagree every iteration)
    auto refresh = [&]() {
        ++refreshes;
        if (!factor()) return false;
        compute_dual();
        repair_dual(true);
        compute_primal();
        fresh = true;
        double obj = 0;
        for (int j = 0; j < N; ++j) obj += c_[j] * x_[j];
        if (best_obj == -kInf || obj > best_obj + 1e-9 * (1 + std::abs(best_obj))) best_obj = obj, progress_at = iterations;
        return true;
    };

    for (;;) {
        if ((iterations & 15) == 0 && elapsed() > time_limit_s) {
            keep_basis();
            return Outcome::TimeLimit;
        }
        if (work_budget().cap && work_budget().used >= work_budget().cap) {  // --work-limit: deterministic stop (denied entries are not counted)
            message = "work limit";
            keep_basis();
            return Outcome::Fallback;
        }
        if (refreshes > 200 + iterations / 30) {  // refactor storm: hand over to the primal long before the stall bound
            message = "dual refactorization storm";
            keep_basis();
            return Outcome::Fallback;
        }
        if (iterations >= kMaxIterations || iterations - progress_at > kStallIterations) {
            message = iterations >= kMaxIterations ? "dual iteration limit" : "dual stalled";
            keep_basis();
            return Outcome::Fallback;
        }
        // Pricing: the most infeasible basic variable relative to its steepest-edge weight.
        int r = -1;
        double best = 0, infeas_r = 0;
        for (int p = 0; p < m_; ++p) {
            int j = head_[p];
            double v = x_[j] < lo_[j] - kPrimalTol ? lo_[j] - x_[j] : x_[j] > up_[j] + kPrimalTol ? x_[j] - up_[j] : 0.0;
            if (v > 0 && v * v > best * weight[p]) r = p, best = v * v / weight[p], infeas_r = v;
        }
        if (r < 0) {
            if (!fresh) {
                if (!refresh()) break;
                continue;
            }
            if (phase1) {
                phase1 = false;
                phase1_iters = iterations;
                best_obj = -kInf, progress_at = iterations;
                lo_ = tlo_, up_ = tup_;
                if (repair_dual(false) > 0) {
                    message = "dual infeasible after dual phase 1";
                    compute_primal();
                    keep_basis();
                    return Outcome::Fallback;
                }
                compute_primal();
                continue;
            }
            keep_basis();
            return Outcome::Optimal;
        }
        int leaving = head_[r];
        bool to_lower = x_[leaving] < lo_[leaving];
        double sigma = to_lower ? -1.0 : 1.0;  // the dual step is t = sigma * theta, theta >= 0

        // Pivot row alpha_rj = rho' a_j for the nonbasics.
        unit.assign(m_, 0.0);
        unit[r] = 1.0;
        btran(unit, rho);
        for (int j : touched) alpha_row[j] = 0, mark[j] = 0;
        touched.clear();
        int nz = 0;
        for (int i = 0; i < m_; ++i) nz += rho[i] != 0;
        if (nz < m_ / 10) {
            for (int i = 0; i < m_; ++i) {
                if (rho[i] == 0) continue;
                for (const Entry& e : rows_[i]) {
                    if (where_[e.index] == kBasic) continue;
                    if (!mark[e.index]) mark[e.index] = 1, touched.push_back(e.index);
                    alpha_row[e.index] += e.value * rho[i];
                }
                int j = n_ + i;
                if (where_[j] != kBasic) {
                    if (!mark[j]) mark[j] = 1, touched.push_back(j);
                    alpha_row[j] -= rho[i];
                }
            }
        } else {
            for (int j = 0; j < N; ++j) {
                if (where_[j] == kBasic) continue;
                double a = dot_col(j, rho);
                if (a != 0) alpha_row[j] = a, mark[j] = 1, touched.push_back(j);
            }
        }

        // Bound-flipping ratio test with Harris tolerances. d_j(theta) = d_j - sigma*theta*alpha_rj.
        cands.clear();
        for (int j : touched) {
            if (lo_[j] == up_[j]) continue;
            double a = sigma * alpha_row[j];
            if (std::abs(a) < kPivotTol) continue;
            bool ok = where_[j] == kZero || (where_[j] == kLower && a > 0) || (where_[j] == kUpper && a < 0);
            if (!ok) continue;
            double ratio = std::max(0.0, d_[j] / a);
            double harris = std::max(0.0, (d_[j] + (a > 0 ? kDualTol : -kDualTol)) / a);
            cands.push_back({j, a, ratio, harris});
        }
        double slope = infeas_r;
        int q = -1;
        size_t bp = 0;  // cands[bp..] are the breakpoints not yet passed (unordered)
        flips.clear();
        while (bp < cands.size()) {
            double hmax = kInf;
            for (size_t k = bp; k < cands.size(); ++k) hmax = std::min(hmax, cands[k].harris);
            size_t end = size_t(std::partition(cands.begin() + bp, cands.end(),
                                               [&](const Cand& c) { return c.ratio <= hmax; }) - cands.begin());
            if (end == bp) end = bp + 1;
            size_t pick = bp;
            double drop = 0;
            for (size_t k = bp; k < end; ++k) {
                if (std::abs(cands[k].a) > std::abs(cands[pick].a)) pick = k;
                drop += (up_[cands[k].j] - lo_[cands[k].j]) * std::abs(cands[k].a);
            }
            if (!(slope - drop > kPrimalTol)) {  // includes an infinite range
                q = cands[pick].j;
                break;
            }
            slope -= drop;
            for (size_t k = bp; k < end; ++k) flips.push_back(cands[k].j);
            bp = end;
        }
        if (q < 0) {
            if (!fresh) {
                if (!refresh()) break;
                continue;
            }
            if (!phase1 && proves_infeasible(r, alpha_row, touched)) {
                keep_basis();
                return Outcome::Infeasible;
            }
            message = phase1 ? "dual phase 1 without an entering variable" : "dual ray not confirmed";
            keep_basis();
            return Outcome::Fallback;
        }

        // Entering column and a stability check against the row computation.
        std::fill(col.begin(), col.end(), 0.0);
        for_col(q, [&](int i, double v) { col[i] = v; });
        ftran(col, alpha);
        double arq = alpha[r];
        if (std::abs(arq - alpha_row[q]) > 1e-7 * (1 + std::abs(arq)) || std::abs(arq) < 1e-11) {
            if (!fresh) {
                if (!refresh()) break;
                continue;
            }
            if (std::abs(arq) < 1e-11) {
                message = "dual pivot too small";
                keep_basis();
                return Outcome::Fallback;
            }
        }
        ++iterations;
        ++work_budget().used;  // work_used counts performed iterations

        // Bound flips: x_B -= B^{-1} sum a_j dx_j.
        if (!flips.empty()) {
            std::fill(delta.begin(), delta.end(), 0.0);
            for (int j : flips) {
                double dx;
                if (where_[j] == kLower) where_[j] = kUpper, dx = up_[j] - lo_[j], x_[j] = up_[j];
                else where_[j] = kLower, dx = lo_[j] - up_[j], x_[j] = lo_[j];
                for_col(j, [&](int i, double v) { delta[i] += v * dx; });
            }
            std::vector<double> dxb;
            ftran(delta, dxb);
            for (int p = 0; p < m_; ++p) x_[head_[p]] -= dxb[p];
        }

        // Dual update: d_j -= t alpha_rj with t = d_q / alpha_rq; the leaving variable gets -t.
        double dq = d_[q];
        if (dq * sigma * alpha_row[q] < 0) {  // slightly infeasible entering reduced cost: shift it to zero
            c_[q] -= dq;
            dq = 0;
            ++shifts;
        }
        double t = dq / arq;
        for (int j : touched)
            if (where_[j] != kBasic) d_[j] -= t * alpha_row[j];
        d_[q] = 0;
        d_[leaving] = -t;

        // Dual steepest-edge weights (Forrest-Goldfarb), with tau = B^{-1} rho.
        double wr = 0;
        for (double v : rho) wr += v * v;
        tau = rho;
        {
            std::vector<double> tmp = tau;
            ftran(tmp, tau);
        }
        for (int p = 0; p < m_; ++p) {
            if (p == r || alpha[p] == 0) continue;
            double kappa = alpha[p] / arq;
            weight[p] = std::max(weight[p] + kappa * (kappa * wr - 2 * tau[p]), 1e-6);
        }
        weight[r] = std::max(wr / (arq * arq), 1e-6);

        // Primal update.
        double rbound = to_lower ? lo_[leaving] : up_[leaving];
        double thp = (x_[leaving] - rbound) / arq;
        for (int p = 0; p < m_; ++p) x_[head_[p]] -= thp * alpha[p];
        x_[q] += thp;
        x_[leaving] = rbound;
        where_[leaving] = to_lower ? kLower : kUpper;
        Eta e{r, arq, {}};
        for (int p = 0; p < m_; ++p)
            if (p != r && alpha[p] != 0) e.col.push_back({p, alpha[p]});
        etas_.push_back(std::move(e));
        head_[r] = q;
        where_[q] = kBasic;
        fresh = false;
        if (etas_.size() >= kRefactorEvery && !refresh()) break;
    }
    keep_basis();
    if (timed_out_) return Outcome::TimeLimit;
    message = "basis could not be refactored";
    return Outcome::Fallback;
}

}  // namespace

Result solve_lp_dual(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                     const std::vector<char>* warm_basis, double time_limit_s) {
    const Stopwatch sw;
    const Model* mp = &model;
    Model neg;
    if (model.maximize) {  // maximise f  ==  minimise -f
        neg = model;
        for (double& c : neg.cost) c = -c;
        neg.obj_const = -neg.obj_const;
        neg.maximize = false;
        mp = &neg;
    }
    Dual dual(*mp, col_lo, col_up);
    Outcome out = dual.run(time_limit_s, warm_basis);
    Result r;
    // The dual engine does not export its ray. Re-solve with primal phase 1 to
    // obtain a checkable proof; exhausted budget returns time_limit, not infeasible.
    if (out == Outcome::TimeLimit) {
        r.status = Status::TimeLimit;
        r.iterations = dual.iterations;
        r.message = dual.note;
        return r;
    }
    // Primal clean-up from the dual basis with the true costs; it certifies the answer.
    double left = time_limit_s - sw();
    r = solve_lp(*mp, col_lo, col_up, dual.basis.empty() ? warm_basis : &dual.basis, left);
    std::string what = "dual " + std::to_string(dual.iterations) + " + primal " + std::to_string(r.iterations) +
                       " iterations" + (dual.message.empty() ? "" : " (" + dual.message + ")") +
                       (dual.note.empty() ? "" : "; " + dual.note) +
                       (dual.start_diag.empty() ? "" : "; " + dual.start_diag + ", phase 1 " + std::to_string(dual.phase1_iters) + " iterations") +
                       "; dual factor " + std::to_string(dual.factor_s) + " s in " + std::to_string(dual.factors) +
                       " (first " + std::to_string(dual.first_factor_s) + " s)";
    r.message = r.message.empty() ? what : r.message + "; " + what;
    r.iterations += dual.iterations;
    if (model.maximize) {
        r.objective = -r.objective;
        r.dual_objective = -r.dual_objective;
        for (double& v : r.row_dual) v = -v;
        for (double& v : r.reduced_cost) v = -v;
    }
    return r;
}
