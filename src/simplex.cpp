// Bounded-variable primal revised simplex with a composite phase 1 (minimise the sum of
// infeasibilities, then the true cost), Devex pricing, a Harris two-pass ratio test,
// bound flips, product-form (eta) basis updates and periodic refactorization.
// Variables 0..n-1 are structural; n+i is the logical of row i with column -e_i, so its
// value equals the activity of row i and its bounds are the row bounds.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>

#include "taral.hpp"

namespace {
constexpr double kPrimalTol = 1e-9;
constexpr double kDualTol = 1e-9;
constexpr double kPivotTol = 1e-9;  // relative to the largest |alpha|
constexpr size_t kRefactorEvery = 100;
constexpr int kPerturbAfter = 3000;    // iterations without objective progress count as a stall
constexpr int kMaxPerturbations = 6;   // stall remedies per solve, so perturb/restore cannot loop forever
constexpr bool kPerturbAtStart = false;
constexpr long kMaxIterations = 50'000'000;

enum Where : char { kBasic, kLower, kUpper, kZero };  // kZero: free nonbasic held at 0

struct Eta {
    int r;
    double pivot;
    std::vector<Entry> col;  // entering column in the old basis, position r excluded
};

class Simplex {
public:
    // clo/cup replace the model's column bounds (solve() passes the model's own).
    Simplex(const Model& md, const std::vector<double>& clo, const std::vector<double>& cup)
        : md_(md), clo_(clo), cup_(cup), m_(int(md.row_names.size())), n_(int(md.col_names.size())) {}
    Result run(double time_limit_s, const std::vector<char>* warm = nullptr);

private:
    const Model& md_;
    const std::vector<double>&clo_, &cup_;
    int m_, n_;
    std::vector<double> lo_, up_, c_, x_;
    std::vector<int> head_;
    std::vector<Where> where_;
    SparseLU lu_;
    std::vector<Eta> etas_;

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
    void set_nonbasic(int j) {
        bool fl = std::isfinite(lo_[j]), fu = std::isfinite(up_[j]);
        if (fl && (!fu || std::abs(x_[j] - lo_[j]) <= std::abs(x_[j] - up_[j]))) where_[j] = kLower, x_[j] = lo_[j];
        else if (fu) where_[j] = kUpper, x_[j] = up_[j];
        else where_[j] = kZero, x_[j] = 0;
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
    // Factor the current basis (repairing rank deficiency with logicals) and recompute x_B.
    bool refactor() {
        for (int attempt = 0; attempt < 5; ++attempt) {
            std::vector<std::vector<Entry>> cols(m_);
            for (int p = 0; p < m_; ++p) for_col(head_[p], [&](int i, double v) { cols[p].push_back({i, v}); });
            etas_.clear();
            if (lu_.factor(m_, cols)) {
                std::vector<double> rhs(m_, 0.0), xb;
                for (int j = 0; j < n_ + m_; ++j)
                    if (where_[j] != kBasic && x_[j] != 0) for_col(j, [&](int i, double v) { rhs[i] -= v * x_[j]; });
                ftran(rhs, xb);
                for (int p = 0; p < m_; ++p) x_[head_[p]] = xb[p];
                return true;
            }
            for (int p : lu_.bad_pos) set_nonbasic(head_[p]);
            for (size_t t = 0; t < lu_.bad_pos.size(); ++t) {
                int j = n_ + lu_.bad_rows[t];
                head_[lu_.bad_pos[t]] = j;
                where_[j] = kBasic;
            }
        }
        return false;
    }
    Result finish(Result r);
};

Result Simplex::finish(Result r) {
    r.x.assign(x_.begin(), x_.begin() + n_);
    r.objective = md_.obj_const;
    for (int j = 0; j < n_; ++j) r.objective += md_.cost[j] * r.x[j];
    // Never claim optimality for a point that fails its own feasibility check.
    std::vector<double> act(m_, 0.0);
    for (int j = 0; j < n_; ++j)
        for (const Entry& e : md_.cols[j]) act[e.index] += e.value * r.x[j];
    double worst = 0;
    for (int i = 0; i < m_; ++i) {
        double scale = 1 + std::abs(std::isfinite(md_.row_lo[i]) ? md_.row_lo[i] : md_.row_up[i]);
        worst = std::max({worst, (md_.row_lo[i] - act[i]) / scale, (act[i] - md_.row_up[i]) / scale});
    }
    for (int j = 0; j < n_; ++j) worst = std::max({worst, clo_[j] - r.x[j], r.x[j] - cup_[j]});
    if (worst > 1e-7) {
        r.status = Status::NumericalFailure;
        r.message = "final point violates constraints by " + std::to_string(worst);
    }
    r.basis.assign(where_.begin(), where_.end());
    return r;
}

Result Simplex::run(double time_limit_s, const std::vector<char>* warm) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    Result res;
    int N = n_ + m_;
    lo_.assign(N, 0), up_.assign(N, 0), c_.assign(N, 0), x_.assign(N, 0);
    where_.assign(N, kBasic), head_.resize(m_);
    for (int j = 0; j < n_; ++j) lo_[j] = clo_[j], up_[j] = cup_[j], c_[j] = md_.cost[j];
    for (int i = 0; i < m_; ++i) lo_[n_ + i] = md_.row_lo[i], up_[n_ + i] = md_.row_up[i], head_[i] = n_ + i;
    for (int j = 0; j < N; ++j)
        if (lo_[j] > up_[j] + kPrimalTol) {
            res.status = Status::Infeasible;
            res.message = "inconsistent bounds";
            return res;
        }
    if (warm && int(warm->size()) == N && std::count(warm->begin(), warm->end(), char(kBasic)) == m_) {
        // Warm start: nonbasics go to their (possibly new) bounds; basics violating a tightened
        // bound are repaired by the composite phase 1.
        for (int j = 0, p = 0; j < N; ++j) {
            where_[j] = Where((*warm)[j]);
            if (where_[j] == kBasic) head_[p++] = j;
            else if (where_[j] == kLower && std::isfinite(lo_[j])) x_[j] = lo_[j];
            else if (where_[j] == kUpper && std::isfinite(up_[j])) x_[j] = up_[j];
            else set_nonbasic(j);
        }
    } else {
        for (int j = 0; j < n_; ++j) set_nonbasic(j);
    }
    if (!refactor()) {
        res.message = "initial basis could not be factored";
        return res;
    }

    bool fresh = true;  // no eta updates since the last refactor
    int stalled = 0;              // iterations since either phase objective last improved
    double best[2] = {kInf, kInf};  // best phase-2 / phase-1 objective seen; flipping phases is no progress
    int level = 0, perturbations = 0;
    // Anti-stalling (idea from the M3 reference engine): widen every non-fixed bound by a tiny
    // deterministic random amount so tied ratios separate; a stall while perturbed retries with
    // a fresh pattern ten times larger. True bounds (level 0) are restored at the perturbed
    // optimum and the solve continues from that basis, so the answer is unperturbed.
    auto apply_bounds = [&](int lvl) {
        std::mt19937_64 rng(12345 + perturbations);
        double scale = lvl ? std::pow(10.0, lvl - 1) : 0.0;
        for (int j = 0; j < N; ++j) {
            double lo = j < n_ ? clo_[j] : md_.row_lo[j - n_];
            double up = j < n_ ? cup_[j] : md_.row_up[j - n_];
            if (lvl && lo != up) {
                double shift = scale * (1e-7 + 1e-6 * double(rng() % 1000000) / 1e6);
                if (std::isfinite(lo)) lo -= shift * (1 + std::abs(lo));
                if (std::isfinite(up)) up += shift * (1 + std::abs(up));
            }
            lo_[j] = lo, up_[j] = up;
            if (where_[j] == kLower) x_[j] = lo;
            else if (where_[j] == kUpper) x_[j] = up;
        }
        level = lvl;
        stalled = 0, best[0] = best[1] = kInf;
        return refactor();  // also recomputes x_B
    };
    if (kPerturbAtStart) {
        ++perturbations;
        if (!apply_bounds(1)) {
            res.message = "perturbed basis could not be factored";
            return res;
        }
    }
    std::vector<double> cb(m_), y, alpha, col(m_), unit, rho;
    std::vector<double> weight(N, 1.0);  // Devex reference weights (approximate edge lengths squared)
    for (;;) {
        if (elapsed() > time_limit_s) {
            res.status = Status::TimeLimit;
            return res;
        }
        if (res.iterations >= kMaxIterations) {
            res.status = Status::IterationLimit;
            return res;
        }
        bool phase1 = false;
        for (int p = 0; p < m_; ++p) {
            int j = head_[p];
            cb[p] = x_[j] < lo_[j] - kPrimalTol ? -1.0 : x_[j] > up_[j] + kPrimalTol ? 1.0 : 0.0;
            phase1 = phase1 || cb[p] != 0;
        }
        if (!phase1)
            for (int p = 0; p < m_; ++p) cb[p] = c_[head_[p]];
        double obj = 0;  // phase objective: sum of infeasibilities, or the true cost
        if (phase1) {
            for (int p = 0; p < m_; ++p) {
                int j = head_[p];
                if (cb[p] < 0) obj += lo_[j] - x_[j];
                else if (cb[p] > 0) obj += x_[j] - up_[j];
            }
        } else {
            for (int j = 0; j < N; ++j) obj += c_[j] * x_[j];
        }
        double& b = best[phase1];
        if (!(b < kInf) || obj < b - 1e-12 * (1 + std::abs(b))) b = obj, stalled = 0;
        else ++stalled;
        if (stalled > kPerturbAfter && perturbations < kMaxPerturbations) {
            ++perturbations;
            if (!apply_bounds(std::min(level + 1, 3))) break;
            fresh = true;
            continue;
        }
        btran(cb, y);

        int q = -1;
        double dq = 0, best_score = 0;
        for (int j = 0; j < N; ++j) {
            if (where_[j] == kBasic || lo_[j] == up_[j]) continue;
            double d = (phase1 ? 0.0 : c_[j]) - dot_col(j, y);
            bool improving = (where_[j] == kLower && d < -kDualTol) || (where_[j] == kUpper && d > kDualTol) ||
                             (where_[j] == kZero && std::abs(d) > kDualTol);
            if (improving && d * d > best_score * weight[j]) q = j, dq = d, best_score = d * d / weight[j];
        }
        if (q < 0) {
            if (!fresh) {
                if (!refactor()) break;
                fresh = true;
                continue;
            }
            if (level > 0) {
                if (!apply_bounds(0)) break;
                continue;
            }
            if (phase1) {
                res.status = Status::Infeasible;
                return res;
            }
            res.status = Status::Optimal;
            return finish(res);
        }

        double dir = dq < 0 ? 1.0 : -1.0;
        std::fill(col.begin(), col.end(), 0.0);
        for_col(q, [&](int i, double v) { col[i] = v; });
        ftran(col, alpha);

        // Harris ratio test. Basic j moves at rate g = -dir*alpha[p] per unit step.
        double amax = 0;
        for (double a : alpha) amax = std::max(amax, std::abs(a));
        double ptol = kPivotTol * std::max(1.0, amax);
        auto blocking = [&](int p, double& g) {  // bound reached by basic at position p, or ±inf
            int j = head_[p];
            g = -dir * alpha[p];
            if (g < 0) return x_[j] > up_[j] + kPrimalTol ? up_[j] : x_[j] >= lo_[j] - kPrimalTol ? lo_[j] : -kInf;
            return x_[j] < lo_[j] - kPrimalTol ? lo_[j] : x_[j] <= up_[j] + kPrimalTol ? up_[j] : kInf;
        };
        double relaxed = kInf;
        for (int p = 0; p < m_; ++p) {
            if (std::abs(alpha[p]) < ptol) continue;
            double g, b = blocking(p, g);
            if (std::isfinite(b)) relaxed = std::min(relaxed, (std::abs(x_[head_[p]] - b) + kPrimalTol) / std::abs(g));
        }
        int r = -1;
        double theta = kInf, rbound = 0;
        for (int p = 0; p < m_; ++p) {
            if (std::abs(alpha[p]) < ptol) continue;
            double g, b = blocking(p, g);
            if (!std::isfinite(b)) continue;
            double ratio = std::max(0.0, (b - x_[head_[p]]) / g);
            if (ratio <= relaxed && (r < 0 || std::abs(alpha[p]) > std::abs(alpha[r]))) r = p, theta = ratio, rbound = b;
        }
        double range = up_[q] - lo_[q];  // inf unless boxed
        if (r < 0 && !std::isfinite(range)) {
            if (!fresh) {
                if (!refactor()) break;
                fresh = true;
                continue;
            }
            if (phase1) {
                res.message = "phase 1 direction without a blocking variable";
                break;
            }
            res.status = Status::Unbounded;
            return res;
        }
        ++res.iterations;
        if (range <= theta) {  // entering variable reaches its opposite bound first
            for (int p = 0; p < m_; ++p) x_[head_[p]] -= dir * range * alpha[p];
            bool to_upper = where_[q] == kLower;
            where_[q] = to_upper ? kUpper : kLower;
            x_[q] = to_upper ? up_[q] : lo_[q];
            continue;
        }
        {  // Devex weight update from the pivot row of the outgoing basis
            unit.assign(m_, 0.0);
            unit[r] = 1.0;
            btran(unit, rho);
            double ar = alpha[r], wq = weight[q];
            for (int j = 0; j < N; ++j) {
                if (where_[j] == kBasic || j == q) continue;
                double a = dot_col(j, rho) / ar;
                if (a != 0) weight[j] = std::max(weight[j], a * a * wq);
            }
            weight[head_[r]] = std::max(wq / (ar * ar), 1.0);
            if (weight[head_[r]] > 1e6) std::fill(weight.begin(), weight.end(), 1.0);  // new reference framework
        }
        for (int p = 0; p < m_; ++p) x_[head_[p]] -= dir * theta * alpha[p];
        x_[q] += dir * theta;
        int leave = head_[r];
        x_[leave] = rbound;
        where_[leave] = rbound == lo_[leave] ? kLower : kUpper;
        Eta e{r, alpha[r], {}};
        for (int p = 0; p < m_; ++p)
            if (p != r && alpha[p] != 0) e.col.push_back({p, alpha[p]});
        etas_.push_back(std::move(e));
        head_[r] = q;
        where_[q] = kBasic;
        fresh = false;
        if (etas_.size() >= kRefactorEvery) {
            if (!refactor()) break;
            fresh = true;
        }
    }
    res.status = Status::NumericalFailure;
    if (res.message.empty()) res.message = "basis could not be refactored";
    return res;
}

}  // namespace

Result solve_lp(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                const std::vector<char>* warm_basis, double time_limit_s) {
    if (!model.maximize) return Simplex(model, col_lo, col_up).run(time_limit_s, warm_basis);
    Model neg = model;  // maximise f  ==  minimise -f
    for (double& c : neg.cost) c = -c;
    neg.obj_const = -neg.obj_const;
    Result r = Simplex(neg, col_lo, col_up).run(time_limit_s, warm_basis);
    r.objective = -r.objective;
    return r;
}

Result solve(const Model& model, double time_limit_s) {
    return solve_lp(model, model.col_lo, model.col_up, nullptr, time_limit_s);
}

const char* status_name(Status s) {
    switch (s) {
        case Status::Optimal: return "optimal";
        case Status::Infeasible: return "infeasible";
        case Status::Unbounded: return "unbounded";
        case Status::TimeLimit: return "time_limit";
        case Status::IterationLimit: return "iteration_limit";
        case Status::NumericalFailure: break;
    }
    return "numerical_failure";
}
