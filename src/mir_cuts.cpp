#include "mir_cuts.hpp"
#include "exact_dyadic.hpp"
#include <algorithm>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <map>

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__)
#error "mir_cuts must not be built with -ffast-math / finite-math / associative / reciprocal math"
#endif
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif
#if defined(__FLT_EVAL_METHOD__) && __FLT_EVAL_METHOD__ != 0
#error "mir_cuts requires FLT_EVAL_METHOD == 0 (no extended-precision intermediates)"
#endif

namespace {
using Clock = std::chrono::steady_clock;
constexpr double kHuge = 1e290;       // anything larger is rejected (overflow guard)
constexpr double kMaxInt = 4503599627370496.0; // 2^52: integer bounds beyond this are not used
constexpr double kMaxMag = 1099511627776.0;    // 2^40: floor/frac operands

struct Ctx { bool bad = false; };
inline double dn(double x) { return std::nextafter(x, -INFINITY); }
inline double upx(double x) { return std::nextafter(x, INFINITY); }
// An fma residual is exactly representable (so "fma(...) == 0" really means exact) only away from underflow. The
// exactness shortcuts below are used only when every operand and the result lie in [2^-300, 2^300]; there the exact
// residual is a multiple of at least 2^-720 and magnitude <= |result|, hence a normal double. Outside that range the
// result is always widened one step outward (never claimed exact).
inline bool safe_range(double v) { double a = std::abs(v); return a >= 0x1p-300 && a <= 0x1p300; }
inline void chk(double v, Ctx& c) { if (!std::isfinite(v) || std::abs(v) > kHuge) c.bad = true; }

// Interval [lo, hi] that encloses the true real value. Each operation is computed in round-to-nearest, so the true
// result of one operation lies within half a spacing of the computed one; when the operation is provably exact (error-free
// transformation or fma) no widening is applied, otherwise the endpoint moves one step outward.
struct Iv { double lo, hi; };
inline Iv iv(double x) { return {x, x}; }
inline Iv neg(Iv a) { return {-a.hi, -a.lo}; }

inline double sum_dir(double a, double b, bool down, Ctx& c) {
    double s = a + b;
    chk(s, c);
    if (c.bad) return s;
    double bb = s - a, err = (a - (s - bb)) + (b - bb);
    if (err == 0.0) return s;
    return down ? dn(s) : upx(s);
}
inline Iv add(Iv a, Iv b, Ctx& c) { return {sum_dir(a.lo, b.lo, true, c), sum_dir(a.hi, b.hi, false, c)}; }
inline Iv sub(Iv a, Iv b, Ctx& c) { return add(a, neg(b), c); }

inline double prod_dir(double x, double y, bool down, Ctx& c) {
    double p = x * y;
    chk(p, c);
    if (c.bad) return p;
    if (x == 0.0 || y == 0.0) return 0.0;
    if (safe_range(x) && safe_range(y) && safe_range(p) && std::fma(x, y, -p) == 0.0) return p;
    return down ? dn(p) : upx(p);
}
inline Iv mul(Iv a, Iv b, Ctx& c) {
    double ls[4], hs[4];
    const double xs[2] = {a.lo, a.hi}, ys[2] = {b.lo, b.hi};
    int k = 0;
    for (double x : xs) for (double y : ys) { ls[k] = prod_dir(x, y, true, c); hs[k] = prod_dir(x, y, false, c); ++k; }
    return {*std::min_element(ls, ls + 4), *std::max_element(hs, hs + 4)};
}
inline double quot_dir(double x, double y, bool down, Ctx& c) {
    double q = x / y;
    chk(q, c);
    if (c.bad) return q;
    if (x == 0.0) return 0.0;
    if (safe_range(x) && safe_range(y) && safe_range(q) && std::fma(q, y, -x) == 0.0) return q;
    return down ? dn(q) : upx(q);
}
inline Iv divi(Iv a, Iv b, Ctx& c) {
    if (!(b.lo > 0 || b.hi < 0)) { c.bad = true; return a; }
    double ls[4], hs[4];
    const double xs[2] = {a.lo, a.hi}, ys[2] = {b.lo, b.hi};
    int k = 0;
    for (double x : xs) for (double y : ys) { ls[k] = quot_dir(x, y, true, c); hs[k] = quot_dir(x, y, false, c); ++k; }
    return {*std::min_element(ls, ls + 4), *std::max_element(hs, hs + 4)};
}
inline double width(Iv a, Ctx& c) { return sum_dir(a.hi, -a.lo, false, c); } // upper bound on hi - lo, >= 0

struct Vub { int y; Iv c0, u; }; // x <= c0 + u*y, real enclosures, y integer with finite integer bounds

struct Prep {
    int m = 0, n = 0;
    std::vector<std::vector<Entry>> rows;
    std::vector<std::vector<Vub>> vub;
    bool timed_out = false;
};

bool int_bounds_ok(double lo, double up) {
    return std::isfinite(lo) && std::isfinite(up) && lo <= up && std::abs(lo) <= kMaxInt && std::abs(up) <= kMaxInt &&
           lo == std::floor(lo) && up == std::floor(up);
}

Prep build(const Model& md, const std::vector<double>& lo, const std::vector<double>& up, int nsrc, Clock::time_point dl) {
    Prep p;
    p.n = int(md.cols.size());
    p.m = nsrc;
    p.rows.assign(p.m, {});
    long seen = 0;
    for (int j = 0; j < p.n; ++j) {
        if ((j & 1023) == 0 && Clock::now() >= dl) { p.timed_out = true; return p; }
        for (const Entry& e : md.cols[j]) {
            if ((++seen & 4095) == 0 && Clock::now() >= dl) { p.timed_out = true; return p; }
            if (e.index < nsrc) p.rows[e.index].push_back({j, e.value});
        }
    }
    p.vub.assign(p.n, {});
    for (int i = 0; i < p.m; ++i) {
        if ((i & 1023) == 0 && Clock::now() >= dl) { p.timed_out = true; return p; }
        if (p.rows[i].size() != 2) continue;
        for (int side = 0; side < 2; ++side) {
            double rhs = side == 0 ? md.row_up[i] : -md.row_lo[i];
            if (!std::isfinite(rhs)) continue;
            double s = side == 0 ? 1.0 : -1.0;
            for (int t = 0; t < 2; ++t) {
                const Entry& ex = p.rows[i][t];
                const Entry& ey = p.rows[i][1 - t];
                int xc = ex.index, yc = ey.index;
                if (md.is_int[xc] || !md.is_int[yc] || !int_bounds_ok(lo[yc], up[yc])) continue;
                double a = s * ex.value, b = s * ey.value;
                if (!(a > 0) || !std::isfinite(a) || !std::isfinite(b)) continue;
                Ctx c;
                Iv c0 = divi(iv(rhs), iv(a), c);
                Iv u = divi(neg(iv(b)), iv(a), c);
                if (c.bad) continue;
                p.vub[xc].push_back({yc, c0, u});
            }
        }
    }
    return p;
}

[[maybe_unused]] double frac_exact(double v) { return v - std::floor(v); } // exact for finite |v| < 2^52

struct SubTerm { int col; Iv c; int kind; Vub v; };
struct IntT { int col; Iv g; int comp; double bnd; double range; };

// Replace each uncertain coefficient by a double so that the new row is implied by the true row:
// int coefficients (variables in [0, range]) take g.lo and their uncertainty is charged to the rhs; continuous s >= 0
// coefficients take c.lo; the rhs becomes the upper endpoint. Returns false (fail closed) on any non-finite value.
bool weaken(const std::vector<Iv>& g, const std::vector<double>& range, const std::vector<Iv>& cc, Iv B, Ctx& c,
            std::vector<double>& gh, std::vector<double>& ch, double& Bh) {
    Iv acc = B;
    gh.resize(g.size());
    for (size_t k = 0; k < g.size(); ++k) {
        gh[k] = g[k].lo;
        double w = width(g[k], c);
        if (w > 0) acc = add(acc, mul(iv(w), iv(range[k]), c), c);
    }
    ch.resize(cc.size());
    for (size_t k = 0; k < cc.size(); ++k) ch[k] = cc[k].lo;
    Bh = acc.hi;
    chk(Bh, c);
    return !c.bad;
}

// Core: the row  sum_j R[k].second * x[R[k].first] <= B0  (coefficients and rhs are enclosures of the true real row, which must
// be a valid inequality over the box [lo, up]; a single source row and an aggregated row both qualify).
bool build_cut_core(const Model& md, const Prep& p, const std::vector<double>& lo, const std::vector<double>& up,
                    const std::vector<double>& x, const std::vector<std::pair<int, Iv>>& R, Iv B0, int row, int side,
                    Clock::time_point dl, MirCut& best, MirStats* st) {
    Ctx c;
    Iv B = B0;
    std::map<int, Iv> icoef;
    std::vector<SubTerm> conts;
    for (const auto& e : R) {
        Iv a = e.second;
        int j = e.first;
        if (md.is_int[j]) {
            if (!int_bounds_ok(lo[j], up[j])) return false;
            auto it = icoef.find(j);
            icoef[j] = add(it == icoef.end() ? iv(0.0) : it->second, a, c);
        } else {
            double dlb = std::isfinite(lo[j]) ? std::max(0.0, x[j] - lo[j]) : kInf;
            double dv = kInf;
            const Vub* bv = nullptr;
            for (const Vub& v : p.vub[j]) {
                double mid = 0.5 * (v.c0.lo + v.c0.hi) + 0.5 * (v.u.lo + v.u.hi) * x[v.y];
                double d = std::max(0.0, mid - x[j]);
                if (d < dv) dv = d, bv = &v;
            }
            if (!std::isfinite(dlb) && !std::isfinite(dv)) return false;
            if (dv < dlb && bv) {
                B = sub(B, mul(a, bv->c0, c), c);
                auto it = icoef.find(bv->y);
                icoef[bv->y] = add(it == icoef.end() ? iv(0.0) : it->second, mul(a, bv->u, c), c);
                conts.push_back({j, neg(a), 1, *bv});
            } else {
                B = sub(B, mul(a, iv(lo[j]), c), c);
                conts.push_back({j, a, 0, Vub{}});
            }
        }
        if (c.bad) { if (st) ++st->cuts_rejected_numerics; return false; }
    }
    std::vector<IntT> ints;
    bool any_frac = false;
    for (auto& kv : icoef) {
        int j = kv.first;
        double bnd, yv;
        int comp = (x[j] - lo[j] <= up[j] - x[j]) ? 0 : 1;
        bnd = comp == 0 ? lo[j] : up[j];
        yv = comp == 0 ? x[j] - lo[j] : up[j] - x[j];
        double range = up[j] - lo[j]; // exact: integers within 2^52
        if (yv > 1e-6 && yv < range - 1e-6) any_frac = true;
        B = sub(B, mul(kv.second, iv(bnd), c), c);
        ints.push_back({j, comp == 0 ? kv.second : neg(kv.second), comp, bnd, range});
    }
    if (c.bad) { if (st) ++st->cuts_rejected_numerics; return false; }
    if (ints.empty() || !any_frac) return false;

    // stage 1 weakening: exact double row  sum gh y' + sum ch s <= Bh,  y' in [0, range], s >= 0
    std::vector<Iv> g1, c1;
    std::vector<double> rg;
    for (const IntT& t : ints) g1.push_back(t.g), rg.push_back(t.range);
    for (const SubTerm& t : conts) c1.push_back(t.c);
    std::vector<double> gh, ch;
    double Bh;
    if (!weaken(g1, rg, c1, B, c, gh, ch, Bh)) { if (st) ++st->cuts_rejected_numerics; return false; }

    // delta candidates
    std::vector<double> deltas;
    for (size_t k = 0; k < ints.size(); ++k) {
        double yv = ints[k].comp == 0 ? x[ints[k].col] - lo[ints[k].col] : up[ints[k].col] - x[ints[k].col];
        if (yv > 1e-6 && yv < ints[k].range - 1e-6 && std::abs(gh[k]) > 1e-6 && std::abs(gh[k]) < 1e9) deltas.push_back(std::abs(gh[k]));
    }
    std::sort(deltas.begin(), deltas.end());
    deltas.erase(std::unique(deltas.begin(), deltas.end()), deltas.end());
    if (deltas.size() > 6) deltas.resize(6);
    std::vector<double> cand;
    for (double d : deltas) for (double div : {1.0, 2.0, 4.0, 8.0}) cand.push_back(d / div);

    bool found = false;
    for (double delta : cand) {
        if (Clock::now() >= dl) { if (st) ++st->timeouts; break; }
        Ctx c2;
        std::vector<Iv> g2(ints.size()), cc2(conts.size());
        for (size_t k = 0; k < ints.size(); ++k) g2[k] = divi(iv(gh[k]), iv(delta), c2);
        for (size_t k = 0; k < conts.size(); ++k) cc2[k] = divi(iv(ch[k]), iv(delta), c2);
        Iv B2 = divi(iv(Bh), iv(delta), c2);
        if (c2.bad) continue;
        std::vector<double> gh2, ch2;
        [[maybe_unused]] double bp;
        double bp0;
        if (!weaken(g2, rg, cc2, B2, c2, gh2, ch2, bp0)) continue;
        if (!(std::abs(bp0) <= kMaxMag)) continue;
        // bp0 - floor(bp0) is generally NOT exactly representable (e.g. bp0 = -0.03 gives 0.97...). Round the fractional part
        // UP to a double f0 (>= the true fraction) and use the weaker row rhs bp = F + f0, which must be exact so that
        // frac(bp) = f0 holds exactly; otherwise reject (fail closed).
        double F = std::floor(bp0);
        Iv tfr = sub(iv(bp0), iv(F), c2);
        double f0 = tfr.hi;
        if (c2.bad || f0 < 0.02 || f0 > 0.98) continue;
        {
            double sbp = F + f0, bb = sbp - F, err = (F - (sbp - bb)) + (f0 - bb);
            if (err != 0.0 || std::floor(sbp) != F || sbp - F != f0) continue;
            bp = sbp;
        }
        Iv omf = sub(iv(1.0), iv(f0), c2); // exact enclosure of 1 - f0
        std::vector<Iv> G(ints.size());
        bool ok = true;
        for (size_t k = 0; k < ints.size(); ++k) {
            double a = gh2[k];
            if (!(std::abs(a) <= kMaxMag)) { ok = false; break; }
            double fl = std::floor(a);
            Iv fa = sub(iv(a), iv(fl), c2); // enclosure of frac(a): the subtraction may be inexact
            Iv term = iv(0.0);
            if (fa.lo > f0) term = divi(sub(fa, iv(f0), c2), omf, c2);
            else if (fa.hi > f0) { Iv u = divi(sub(fa, iv(f0), c2), omf, c2); term = {0.0, std::max(0.0, u.hi)}; }
            G[k] = add(iv(fl), term, c2);
        }
        std::vector<Iv> H(conts.size());
        std::vector<char> keep(conts.size(), 0);
        for (size_t k = 0; k < conts.size() && ok; ++k) {
            if (ch2[k] < 0) { H[k] = divi(iv(ch2[k]), omf, c2); keep[k] = 1; }
        }
        if (!ok || c2.bad) continue;
        // stage 3 weakening of the cut coefficients (y' in [0, range], s >= 0)
        std::vector<Iv> Hk;
        std::vector<size_t> hidx;
        for (size_t k = 0; k < conts.size(); ++k) if (keep[k]) Hk.push_back(H[k]), hidx.push_back(k);
        std::vector<double> Gh, Hh;
        double Fh;
        if (!weaken(G, rg, Hk, iv(F), c2, Gh, Hh, Fh)) continue;
        // stage 4: back-transform to x space with enclosures
        std::map<int, Iv> coef;
        auto addc = [&](int j, Iv v) { auto it = coef.find(j); coef[j] = add(it == coef.end() ? iv(0.0) : it->second, v, c2); };
        Iv rhs = iv(Fh);
        for (size_t k = 0; k < ints.size(); ++k) {
            const IntT& t = ints[k];
            if (t.comp == 0) { addc(t.col, iv(Gh[k])); rhs = add(rhs, mul(iv(Gh[k]), iv(t.bnd), c2), c2); }
            else { addc(t.col, iv(-Gh[k])); rhs = sub(rhs, mul(iv(Gh[k]), iv(t.bnd), c2), c2); }
        }
        bool used_vub = false;
        for (size_t q = 0; q < hidx.size(); ++q) {
            const SubTerm& t = conts[hidx[q]];
            double h = Hh[q];
            if (t.kind == 0) { addc(t.col, iv(h)); rhs = add(rhs, mul(iv(h), iv(lo[t.col]), c2), c2); }
            else {
                addc(t.col, iv(-h));
                addc(t.v.y, mul(iv(h), t.v.u, c2));
                rhs = sub(rhs, mul(iv(h), t.v.c0, c2), c2);
                used_vub = true;
            }
        }
        if (c2.bad) { if (st) ++st->cuts_rejected_numerics; continue; }
        // final: pick coefficient doubles, charge the coefficient uncertainty over the magnitude domain of x
        MirCut cut;
        double mx = 0, mn = kInf, act = 0, nrm = 0;
        bool reject = false;
        for (auto& kv : coef) {
            int j = kv.first;
            double v = kv.second.lo;
            double w = width(kv.second, c2);
            if (w > 0) {
                double M = std::max(std::abs(lo[j]), std::abs(up[j]));
                if (!std::isfinite(M)) { reject = true; if (st) ++st->cuts_rejected_domain; break; }
                rhs = add(rhs, mul(iv(w), iv(M), c2), c2);
            }
            if (v == 0.0) continue;
            mx = std::max(mx, std::abs(v)), mn = std::min(mn, std::abs(v));
            act += v * x[j], nrm += v * v;
            cut.idx.push_back(j), cut.val.push_back(v);
        }
        if (reject || c2.bad || cut.idx.empty() || mx / mn > 1e6 || cut.idx.size() > 500) { if (!reject && st) ++st->cuts_rejected_numerics; continue; }
        double r = rhs.hi;
        chk(r, c2);
        if (c2.bad) continue;
        double eff = (act - r) / std::sqrt(nrm);
        if (!(eff > 0)) continue;
        cut.rhs = r, cut.efficacy = eff, cut.source_row = row, cut.side = side;
        if (!found || eff > best.efficacy) { best = cut; found = true; if (st) st->vub_used += used_vub ? 0 : 0; best.side = side + (used_vub ? 2 : 0); }
    }
    return found;
}

// Single source row (unchanged behaviour: the same enclosures as before the refactor).
bool build_cut(const Model& md, const Prep& p, const std::vector<double>& lo, const std::vector<double>& up,
               const std::vector<double>& x, int row, int side, Clock::time_point dl, MirCut& best, MirStats* st) {
    const auto& RR = p.rows[row];
    double s = side == 0 ? 1.0 : -1.0;
    double beta0 = side == 0 ? md.row_up[row] : -md.row_lo[row];
    if (!std::isfinite(beta0)) return false;
    if (RR.size() < 2 || RR.size() > 400) { if (st) ++st->rows_skipped; return false; }
    std::vector<std::pair<int, Iv>> R;
    R.reserve(RR.size());
    for (const Entry& e : RR) R.push_back({e.index, iv(s * e.value)});
    return build_cut_core(md, p, lo, up, x, R, iv(beta0), row, side, dl, best, st);
}

// R14 (opt-in, max_agg > 0): Marchand-Wolsey path aggregation. Starting from a tight source row, repeatedly eliminate a continuous
// variable that is strictly inside its bounds by adding a nonnegative multiple lam of another tight source row (lam is any
// positive double: the real combination is a valid inequality whatever lam is; its coefficients and rhs are enclosed by
// interval arithmetic). The eliminated variable's residual coefficient (an enclosure of ~0) is folded into the rhs through a
// rigorous bound on |x_j| (finite lower bound and a finite simple or VUB-derived upper bound required), so no tiny noise term
// is kept. After each step the core c-MIR runs on the aggregated row. Selection uses double approximations only.
bool agg_cut(const Model& md, const Prep& p, const std::vector<double>& lo, const std::vector<double>& up, const std::vector<double>& x,
             const std::vector<std::vector<int>>& colrows, const std::vector<double>& act, int row, int side, int max_agg,
             Clock::time_point dl, MirCut& best, MirStats* st) {
    double s0 = side == 0 ? 1.0 : -1.0;
    double b0 = side == 0 ? md.row_up[row] : -md.row_lo[row];
    if (!std::isfinite(b0)) return false;
    Ctx c;
    std::map<int, Iv> A;
    for (const Entry& e : p.rows[row]) A[e.index] = iv(s0 * e.value);
    Iv B = iv(b0);
    std::vector<int> used{row};
    // rigorous upper bound of |x_j| (needs finite lower bound and a finite simple or VUB-derived upper bound)
    auto abs_bound = [&](int j, double& M) {
        if (!std::isfinite(lo[j])) return false;
        double ubx = std::isfinite(up[j]) ? up[j] : kInf;
        for (const Vub& v : p.vub[j]) {
            if (!int_bounds_ok(lo[v.y], up[v.y])) continue;
            Ctx cc;
            double hi = -kInf;
            for (double ye : {lo[v.y], up[v.y]}) hi = std::max(hi, add(v.c0, mul(v.u, iv(ye), cc), cc).hi);
            if (!cc.bad) ubx = std::min(ubx, hi);
        }
        if (!std::isfinite(ubx)) return false;
        M = std::max(std::abs(lo[j]), std::abs(ubx));
        return std::isfinite(M) && M < kHuge;
    };
    for (int step = 1; step <= max_agg; ++step) {
        if (Clock::now() >= dl) return false;
        struct Cand { double d; int j; };
        std::vector<Cand> cand;
        for (const auto& kv : A) {
            int j = kv.first;
            if (md.is_int[j]) continue;
            double cm = 0.5 * (kv.second.lo + kv.second.hi);
            if (std::abs(cm) < 1e-12) continue;
            double d = std::max(0.0, x[j] - lo[j]);
            if (std::isfinite(up[j])) d = std::min(d, std::max(0.0, up[j] - x[j]));
            for (const Vub& v : p.vub[j]) d = std::min(d, std::max(0.0, 0.5 * (v.c0.lo + v.c0.hi) + 0.5 * (v.u.lo + v.u.hi) * x[v.y] - x[j]));
            if (d > 1e-6) cand.push_back({d, j});
        }
        if (cand.empty()) return false;
        std::sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b) { return a.d > b.d || (a.d == b.d && a.j < b.j); });
        bool done = false;
        for (size_t t = 0; t < cand.size() && t < 3 && !done; ++t) {
            int j = cand[t].j;
            double cm = 0.5 * (A[j].lo + A[j].hi);
            int bi = -1, bs = 0;
            double blam = 0;
            size_t bsz = 0;
            for (int i2 : colrows[j]) {
                if (std::find(used.begin(), used.end(), i2) != used.end()) continue;
                for (int sd = 0; sd < 2; ++sd) {
                    double bd = sd == 0 ? md.row_up[i2] : -md.row_lo[i2];
                    if (!std::isfinite(bd)) continue;
                    double s2 = sd == 0 ? 1.0 : -1.0;
                    if (bd - s2 * act[i2] > 1e-6 * (1 + std::abs(bd))) continue; // row must be tight at x
                    double aj = 0;
                    for (const Entry& e : p.rows[i2]) if (e.index == j) aj = s2 * e.value;
                    if (aj == 0.0) continue;
                    double lam = -cm / aj;
                    if (!(lam > 0) || !std::isfinite(lam) || lam > 1e6 || lam < 1e-6) continue;
                    if (bi < 0 || p.rows[i2].size() < bsz) bi = i2, bs = sd, blam = lam, bsz = p.rows[i2].size();
                }
            }
            if (bi < 0) continue;
            double s2 = bs == 0 ? 1.0 : -1.0, bd = bs == 0 ? md.row_up[bi] : -md.row_lo[bi];
            std::map<int, Iv> A2;
            Iv B2{0, 0};
            bool got = false;
            double lamv[5] = {blam, std::nextafter(blam, INFINITY), std::nextafter(blam, 0.0), std::nextafter(std::nextafter(blam, INFINITY), INFINITY), std::nextafter(std::nextafter(blam, 0.0), 0.0)};
            for (double lv : lamv) {
                Ctx c3;
                Iv lam = iv(lv);
                std::map<int, Iv> At = A;
                for (const Entry& e : p.rows[bi]) {
                    auto it = At.find(e.index);
                    At[e.index] = add(it == At.end() ? iv(0.0) : it->second, mul(lam, iv(s2 * e.value), c3), c3);
                }
                Iv Bt = add(B, mul(lam, iv(bd), c3), c3);
                Iv rj = At[j];
                double M;
                if (rj.lo == 0.0 && rj.hi == 0.0) { /* exact cancellation */ }
                else if (std::isfinite(lo[j]) && lo[j] >= 0.0 && rj.lo >= 0.0) { /* -r*x <= 0 for r,x >= 0 */ }
                else if (abs_bound(j, M)) { double w = std::max(std::abs(rj.lo), std::abs(rj.hi)); Bt = add(Bt, mul(iv(w), iv(M), c3), c3); }
                else continue;
                if (c3.bad) continue;
                At.erase(j);
                A2 = std::move(At), B2 = Bt, got = true;
                break;
            }
            if (!got) continue;
            A = std::move(A2);
            B = B2;
            used.push_back(bi);
            done = true;
        }
        if (!done) return false;
        std::vector<std::pair<int, Iv>> R;
        for (const auto& kv : A) if (!(kv.second.lo == 0.0 && kv.second.hi == 0.0)) R.push_back({kv.first, kv.second});
        if (R.size() < 2 || R.size() > 400) return false;
        if (st) ++st->agg_tried;
        MirCut cc;
        if (build_cut_core(md, p, lo, up, x, R, B, row, side, dl, cc, st)) { best = cc; return true; }
    }
    return false;
}
} // namespace

// a*b+c in one expression: a compiler allowed to contract (-ffp-contract=fast with FMA hardware) returns the fused value.
__attribute__((noinline)) static double contract_probe(double a, double b, double c) { return a * b + c; }

bool fp_assumptions_ok() {
    if (std::fegetround() != FE_TONEAREST) return false;
    volatile double one = 1.0, h = 0x1p-53, t = 0x1.8p-53;
    volatile double a = one + h;       // tie: rounds to even (1.0)
    volatile double b = one + t;       // above the tie: rounds up to 1 + 2^-52
    volatile double big = 9007199254740993.0; // 2^53 + 1 is not representable: rounds to even 2^53
    if (!(a == 1.0 && b == 1.0 + 0x1p-52 && big == 9007199254740992.0 && sizeof(double) == 8 && DBL_MANT_DIG == 53)) return false;
#if defined(__x86_64__) || defined(__i386__)
    unsigned csr = _mm_getcsr();
    if (csr & 0x8040) return false; // flush-to-zero (0x8000) or denormals-are-zero (0x0040)
#endif
    // gradual underflow must really be in force (also covers non-x86): subnormal operands and subnormal results survive
    volatile double tiny = 0x1p-1074, edge = 0x1p-1022, half = 0.5;
    volatile double s1 = tiny + tiny;          // 2^-1073, needs DAZ off
    volatile double s2 = edge * half;          // 2^-1023 subnormal result, needs FTZ off
    volatile double s3 = (edge + tiny) - edge; // 2^-1074: gradual underflow keeps the difference
    if (!(s1 == 0x1p-1073 && s2 == 0x1p-1023 && s3 == 0x1p-1074)) return false;
    // contraction probe: (1+2^-30)(1-2^-30) - 1 is 0 unfused and -2^-60 when fused
    volatile double pa = 1.0 + 0x1p-30, pb = 1.0 - 0x1p-30, pc = -1.0;
    if (contract_probe(pa, pb, pc) != 0.0) return false;
    return true;
}

std::vector<MirCut> separate_mir(const Model& md, const std::vector<double>& lo, const std::vector<double>& up,
                                 const std::vector<double>& x, int nsrc_rows, int max_cuts, double min_eff,
                                 Clock::time_point deadline, MirStats* stats, int max_agg) {
    std::vector<MirCut> out;
    if (!fp_assumptions_ok()) { if (stats) stats->fp_ok = false; return out; }
    if (nsrc_rows > int(md.row_lo.size())) nsrc_rows = int(md.row_lo.size());
    Prep p = build(md, lo, up, nsrc_rows, deadline);
    if (p.timed_out) { if (stats) ++stats->timeouts; return out; }
    std::vector<MirCut> all;
    std::vector<std::vector<int>> colrows;
    std::vector<double> act;
    if (max_agg > 0) {
        colrows.assign(p.n, {});
        act.assign(p.m, 0.0);
        for (int i = 0; i < p.m; ++i)
            for (const Entry& e : p.rows[i]) { colrows[e.index].push_back(i); act[i] += e.value * x[e.index]; }
    }
    for (int i = 0; i < p.m; ++i) {
        if (Clock::now() >= deadline) { if (stats) ++stats->timeouts; break; }
        for (int side = 0; side < 2; ++side) {
            double b = side == 0 ? md.row_up[i] : -md.row_lo[i];
            if (!std::isfinite(b)) continue;
            if (stats) ++stats->rows_tried;
            MirCut c;
            if (build_cut(md, p, lo, up, x, i, side, deadline, c, stats) && c.efficacy >= min_eff) { all.push_back(std::move(c)); continue; }
            if (max_agg > 0 && b - (side == 0 ? 1.0 : -1.0) * act[i] <= 1e-6 * (1 + std::abs(b)) && p.rows[i].size() >= 2) {
                MirCut ca;
                if (agg_cut(md, p, lo, up, x, colrows, act, i, side, max_agg, deadline, ca, stats) && ca.efficacy >= min_eff) {
                    if (stats) ++stats->agg_cuts;
                    all.push_back(std::move(ca));
                }
            }
        }
    }
    std::stable_sort(all.begin(), all.end(), [](const MirCut& a, const MirCut& b) { return a.efficacy > b.efficacy; });
    std::vector<double> norms;
    for (MirCut& c : all) {
        if (int(out.size()) >= max_cuts) break;
        if (Clock::now() >= deadline) { if (stats) ++stats->timeouts; break; }
        bool par = false;
        double nc = 0;
        for (double v : c.val) nc += v * v;
        for (size_t oi = 0; oi < out.size() && !par; ++oi) {
            const MirCut& o = out[oi];
            double dot = 0;
            size_t a = 0, b = 0;
            while (a < c.idx.size() && b < o.idx.size()) {
                if (c.idx[a] == o.idx[b]) dot += c.val[a++] * o.val[b++];
                else if (c.idx[a] < o.idx[b]) ++a;
                else ++b;
            }
            if (std::abs(dot) / std::sqrt(nc * norms[oi]) > 0.98) par = true;
        }
        if (!par) {
            norms.push_back(nc);
            if (stats && c.side >= 2) ++stats->vub_used;
            c.side &= 1;
            out.push_back(std::move(c));
        }
    }
    if (stats) stats->cuts_found += long(out.size());
    return out;
}

CutCheck check_cut_exact(const MirCut& c, const std::vector<double>& x, double tol) {
    using namespace exact_dyadic;
    Wide act, w;
    for (size_t k = 0; k < c.idx.size(); ++k) {
        if (c.idx[k] < 0 || size_t(c.idx[k]) >= x.size() || !std::isfinite(x[c.idx[k]]) || !std::isfinite(c.val[k])) return CutCheck::Inconclusive;
        w = Wide();
        if (!term(w, c.val[k], x[c.idx[k]]) || !add(act, w)) return CutCheck::Inconclusive;
    }
    double slack = c.rhs + tol * (1 + std::abs(c.rhs));
    Wide bound;
    if (!std::isfinite(slack) || !term(bound, slack, 1.0)) return CutCheck::Inconclusive;
    return compare(act, bound) <= 0 ? CutCheck::Holds : CutCheck::Violated;
}
