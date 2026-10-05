// Original-coordinate LP proof checker. No presolve data or basis flags are trusted.
#include <algorithm>
#include <cmath>
#include <limits>
#include "taral.hpp"

static bool verify_impl(const Model& md, const std::vector<double>& lo,
                        const std::vector<double>& up, Result& r) {
    constexpr double tol = 1e-8;
    const size_t n = md.cols.size(), m = md.row_names.size();
    auto vector_ok = [](const std::vector<double>& v, size_t size, bool nonnegative) {
        if (v.size() != size) return false;
        for (double a : v) if (!std::isfinite(a) || (nonnegative && a < 0)) return false;
        return true;
    };
    bool ok = lo.size() == n && up.size() == n;
    double residual = 0, margin = 0;
    ok = ok && md.row_lo.size() == m && md.row_up.size() == m && md.cost.size() == n;
    for (size_t j = 0; ok && j < n; ++j) {
        ok = std::isfinite(md.cost[j]) && !std::isnan(lo[j]) && !std::isnan(up[j]) &&
             lo[j] != kInf && up[j] != -kInf;
        for (const Entry& e : md.cols[j])
            ok = ok && e.index >= 0 && size_t(e.index) < m && std::isfinite(e.value);
    }
    for (size_t i = 0; ok && i < m; ++i)
        ok = !std::isnan(md.row_lo[i]) && !std::isnan(md.row_up[i]) &&
             md.row_lo[i] != kInf && md.row_up[i] != -kInf;
    if (r.status == Status::Infeasible && ok) {
        ok = vector_ok(r.farkas_row_lower, m, true) && vector_ok(r.farkas_row_upper, m, true) &&
             vector_ok(r.farkas_col_lower, n, true) && vector_ok(r.farkas_col_upper, n, true);
        if (ok) {
            std::vector<long double> z(n);
            // Aggregate duplicate (row,column) entries before any bound inference.
            std::vector<std::vector<Entry>> agg(n);
            std::vector<int> row_nnz(m, 0);
            for (size_t j = 0; j < n; ++j) {
                std::vector<Entry> sorted = md.cols[j];
                std::stable_sort(sorted.begin(), sorted.end(),
                                 [](const Entry& a, const Entry& b) { return a.index < b.index; });
                for (size_t k = 0; k < sorted.size();) {
                    size_t q = k; long double sum = 0;
                    for (; q < sorted.size() && sorted[q].index == sorted[k].index; ++q) sum += sorted[q].value;
                    if (sum != 0) { agg[j].push_back({sorted[k].index, (double)sum}); ++row_nnz[sorted[k].index]; }
                    k = q;
                }
            }
            // Implied column bounds from row activities (computed lazily, only when a
            // column with a roundoff stationarity residual has no finite bound). Every
            // feasible point of the node satisfies them, so z*implied_bound is a valid
            // compensation for ALL feasible x, not an assumption about x.
            std::vector<double> ilo, iup;
            bool implied_ready = false;
            auto build_implied = [&]() {
                implied_ready = true;
                ilo = lo; iup = up;
                std::vector<std::vector<Entry>> rows(m);  // Entry.index = column here
                for (size_t j = 0; j < n; ++j)
                    for (const Entry& e : agg[j]) rows[e.index].push_back({(int)j, e.value});
                for (int pass = 0; pass < 8; ++pass) {
                    bool changed = false;
                    for (size_t i = 0; i < m; ++i) {
                        long double mn = 0, mx = 0; int mn_inf = 0, mx_inf = 0;
                        for (const Entry& e : rows[i]) {
                            double l = ilo[e.index], u = iup[e.index];
                            double lowc = e.value > 0 ? l : u, highc = e.value > 0 ? u : l;
                            if (std::isfinite(lowc)) mn += (long double)e.value*lowc; else ++mn_inf;
                            if (std::isfinite(highc)) mx += (long double)e.value*highc; else ++mx_inf;
                        }
                        for (const Entry& e : rows[i]) {
                            double l = ilo[e.index], u = iup[e.index], a = e.value;
                            double lowc = a > 0 ? l : u, highc = a > 0 ? u : l;
                            bool max_rest = mx_inf == 0 || (mx_inf == 1 && !std::isfinite(highc));
                            bool min_rest = mn_inf == 0 || (mn_inf == 1 && !std::isfinite(lowc));
                            long double maxrest = mx - (std::isfinite(highc) ? (long double)a*highc : 0);
                            long double minrest = mn - (std::isfinite(lowc) ? (long double)a*lowc : 0);
                            double nl = -kInf, nu = kInf;  // new bounds on x_j
                            if (std::isfinite(md.row_lo[i]) && max_rest) {
                                double v = (double)((md.row_lo[i] - maxrest)/a);
                                if (a > 0) nl = v; else nu = v;
                            }
                            if (std::isfinite(md.row_up[i]) && min_rest) {
                                double v = (double)((md.row_up[i] - minrest)/a);
                                if (a > 0) nu = std::min(nu, v); else nl = std::max(nl, v);
                            }
                            size_t j = e.index;
                            if (std::isfinite(nl) && nl > ilo[j]) { ilo[j] = nl; changed = true; }
                            if (std::isfinite(nu) && nu < iup[j]) { iup[j] = nu; changed = true; }
                        }
                    }
                    if (!changed) break;
                }
            };
            long double value = 0, scale = 0, norm = 0;
            auto term = [&](double a, double b, int sign) {
                norm += a;
                if (a == 0) return;
                if (!std::isfinite(b)) { ok = false; return; }
                value += sign * (long double)a*b; scale += std::abs((long double)a*b);
            };
            for (size_t i = 0; i < m; ++i) {
                term(r.farkas_row_lower[i], md.row_lo[i], 1);
                term(r.farkas_row_upper[i], md.row_up[i], -1);
            }
            for (size_t j = 0; j < n; ++j) {
                z[j] = -r.farkas_col_lower[j] + r.farkas_col_upper[j];
                for (const Entry& e : md.cols[j])
                    z[j] += (long double)e.value * (r.farkas_row_upper[e.index]-r.farkas_row_lower[e.index]);
                term(r.farkas_col_lower[j], lo[j], 1); term(r.farkas_col_upper[j], up[j], -1);
                // Residual stationarity cannot be ignored on an unbounded column:
                // compensate it using a finite minimizing bound, or reject.
                if (z[j] != 0) {
                    double bound = z[j] > 0 ? lo[j] : up[j];
                    if (!std::isfinite(bound)) {
                        // Exact repair: a residual z>0 can be removed by lowering the column's
                        // upper multiplier by z (z<0: lowering the lower multiplier by |z|). The
                        // result is another valid nonnegative multiplier set with z=0, valid only
                        // when that multiplier covers |z| and the matching bound is finite. Its
                        // objective term changes by z*up (z>0) or z*lo (z<0). No assumption on x
                        // is made, so this stays sound for unbounded columns.
                        double mult = z[j] > 0 ? r.farkas_col_upper[j] : r.farkas_col_lower[j];
                        double other = z[j] > 0 ? up[j] : lo[j];
                        if (std::isfinite(other) && (long double)mult >= std::abs(z[j]))
                            bound = other;
                    }
                    if (!std::isfinite(bound)) {
                        // A one-sided column can still have an ORIGINAL singleton
                        // row bound. Use that logically implied finite bound to
                        // control stationarity roundoff; never simply ignore it.
                        // Use the AGGREGATED (row,column) coefficient, never a
                        // single duplicate entry, and only rows whose aggregated
                        // matrix row has exactly one nonzero column.
                        for (const Entry& e : agg[j]) {
                            if (row_nnz[e.index] != 1) continue;
                            double endpoint = (z[j] > 0) == (e.value > 0) ? md.row_lo[e.index] : md.row_up[e.index];
                            double candidate = endpoint/e.value;
                            if (std::isfinite(candidate)) {
                                if (!std::isfinite(bound)) bound = candidate;
                                else bound = z[j] > 0 ? std::max(bound,candidate) : std::min(bound,candidate);
                            }
                        }
                    }
                    if (!std::isfinite(bound)) {
                        if (!implied_ready) build_implied();
                        double cand = z[j] > 0 ? ilo[j] : iup[j];
                        if (std::isfinite(cand)) bound = cand;
                    }
                    if (!std::isfinite(bound)) ok = false;
                    else { value += z[j]*bound; scale += std::abs(z[j]*bound); }
                }
                residual = std::max(residual, (double)std::abs(z[j]));
            }
            margin = (double)(value/(1+scale));
            ok = ok && std::abs(norm-1) <= tol && residual <= tol && margin > tol;
            // Nonzero stationarity can escape along a free column. Demand near
            // machine precision cancellation, not merely a loose feasibility test.
            ok = ok && residual <= 1e-12;
        }
    } else if (r.status == Status::Unbounded && ok) {
        ok = vector_ok(r.x, n, false) && vector_ok(r.ray, n, false);
        if (ok) {
            std::vector<long double> act(m, 0), direction(m, 0);
            long double slope = 0, cnorm = 0;
            double norm = 0;
            auto check = [&](long double x, long double d, double l, double u) {
                if (std::isfinite(l)) {
                    residual = std::max(residual, (double)std::max((long double)0, (l-x)/(1+std::abs(l))));
                    residual = std::max(residual, (double)std::max((long double)0, -d));
                }
                if (std::isfinite(u)) {
                    residual = std::max(residual, (double)std::max((long double)0, (x-u)/(1+std::abs(u))));
                    residual = std::max(residual, (double)std::max((long double)0, d));
                }
            };
            for (size_t j = 0; j < n; ++j) {
                norm = std::max(norm, std::abs(r.ray[j]));
                slope += (long double)md.cost[j]*r.ray[j]; cnorm += std::abs(md.cost[j]);
                check(r.x[j], r.ray[j], lo[j], up[j]);
                for (const Entry& e : md.cols[j]) {
                    act[e.index] += (long double)e.value*r.x[j];
                    direction[e.index] += (long double)e.value*r.ray[j];
                }
            }
            for (size_t i = 0; i < m; ++i) check(act[i], direction[i], md.row_lo[i], md.row_up[i]);
            margin = (double)((md.maximize ? slope : -slope)/(1+cnorm));
            ok = residual <= tol && std::abs(norm-1) <= tol && margin > tol;
        }
    } else ok = false;
    ok = ok && std::isfinite(residual) && std::isfinite(margin);
    r.certificate_verified = ok;
    r.certificate_residual = std::isfinite(residual) ? residual : 0;
    r.certificate_margin = std::isfinite(margin) ? margin : 0;
    if (!ok) { r.status = Status::NumericalFailure; r.message = "nonoptimal original-model proof failed or unavailable"; }
    return ok;
}

// ---- Strict-contract Farkas multiplier repair (emission side) ----
// The independent strict checker (benchmarks/nonoptimal_certificate_check.py) sums every
// stationarity coordinate exactly (fsum of rounded products) and rejects ANY nonzero value on a
// column whose needed-side bound is not finite (after singleton-row rescue). Roundoff of 1e-19
// is enough to reject. Here the emitted multipliers are repaired by bounded local moves that keep
// every multiplier nonnegative and only touch finite-bound support, then re-tested with a port of
// the checker's own arithmetic. A repair is adopted ONLY when that port accepts it AND the
// original verifier still accepts it; otherwise the original multipliers are kept unchanged.
namespace {
double exact_sum(std::vector<double>& v) {  // port of CPython math.fsum (exactly rounded)
    std::vector<double> p;
    for (double x : v) {
        size_t i = 0;
        for (double y : p) {
            if (std::abs(x) < std::abs(y)) std::swap(x, y);
            double hi = x + y, lo = y - (hi - x);
            if (lo != 0) p[i++] = lo;
            x = hi;
        }
        p.resize(i); p.push_back(x);
    }
    size_t n = p.size(); double hi = 0;
    if (n > 0) {
        hi = p[--n]; double lo = 0;
        while (n > 0) {
            double x = hi, y = p[--n]; hi = x + y; double yr = hi - x; lo = y - yr;
            if (lo != 0) break;
        }
        if (n > 0 && ((lo < 0 && p[n-1] < 0) || (lo > 0 && p[n-1] > 0))) {
            double y = lo * 2, x = hi + y, yr = x - hi;
            if (y == yr) hi = x;
        }
    }
    return hi;
}
struct StrictCtx {
    const Model& md; const std::vector<double>& lo; const std::vector<double>& up;
    size_t n, m; std::vector<std::vector<Entry>> agg; std::vector<int> row_nnz;
    // pyagg[j]: per-column (row,value) with duplicate (row,column) entries summed sequentially in
    // double in file order, exactly like the strict checker's parser (col[r] = col.get(r,0.0)+v).
    std::vector<std::vector<Entry>> pyagg;
    long long work = 0;
    StrictCtx(const Model& md_, const std::vector<double>& l, const std::vector<double>& u)
        : md(md_), lo(l), up(u), n(md_.cols.size()), m(md_.row_names.size()), agg(n), row_nnz(m, 0), pyagg(n) {
        for (size_t j = 0; j < n; ++j) {
            // first-appearance order of rows, sequential double sums (stable on file order)
            std::vector<Entry> s;
            for (const Entry& e : md.cols[j]) {
                bool found = false;
                for (Entry& q : s) if (q.index == e.index) { q.value += e.value; found = true; break; }
                if (!found) s.push_back(e);
            }
            for (const Entry& q : s) { pyagg[j].push_back(q); if (q.value != 0) ++row_nnz[q.index]; }
        }
    }
    // checker-exact stationarity of column j for multiplier set (rl,ru,cl,cu); NaN on overflow
    double zcol(size_t j, const std::vector<double>& rl, const std::vector<double>& ru,
                const std::vector<double>& cl, const std::vector<double>& cu) {
        std::vector<double> t; t.reserve(pyagg[j].size() + 2);
        t.push_back(cu[j]); t.push_back(-cl[j]);
        for (const Entry& e : pyagg[j]) {
            double yv = ru[e.index] - rl[e.index], pr = e.value * yv;
            if (!std::isfinite(yv) || !std::isfinite(pr)) return std::numeric_limits<double>::quiet_NaN();
            t.push_back(pr);
        }
        work += (long long)t.size();
        return exact_sum(t);
    }
    // needed-side bound after the checker's singleton-row rescue; NaN-free, may be infinite
    double need_bound(size_t j, double z) const {
        double bound = z > 0 ? lo[j] : up[j];
        if (std::isfinite(bound)) return bound;
        for (const Entry& e : pyagg[j]) {
            if (e.value == 0 || row_nnz[e.index] != 1) continue;
            double endpoint = (z > 0) == (e.value > 0) ? md.row_lo[e.index] : md.row_up[e.index];
            double cand = endpoint / e.value;
            if (std::isfinite(cand)) bound = !std::isfinite(bound) ? cand : (z > 0 ? std::max(bound, cand) : std::min(bound, cand));
        }
        return bound;
    }
};
}  // namespace

static bool strict_farkas_ok(StrictCtx& C, const Result& r) {
    const size_t n = C.n, m = C.m;
    std::vector<double> norm_t(r.farkas_row_lower); norm_t.insert(norm_t.end(), r.farkas_row_upper.begin(), r.farkas_row_upper.end());
    norm_t.insert(norm_t.end(), r.farkas_col_lower.begin(), r.farkas_col_lower.end());
    norm_t.insert(norm_t.end(), r.farkas_col_upper.begin(), r.farkas_col_upper.end());
    double norm = exact_sum(norm_t); std::vector<double> terms;
    auto term = [&](double a, double b, int sgn) { if (a) { if (!std::isfinite(b)) return false; double t = sgn * a * b; if (!std::isfinite(t)) return false; terms.push_back(t); } return true; };
    for (size_t i = 0; i < m; ++i)
        if (!term(r.farkas_row_lower[i], C.md.row_lo[i], 1) || !term(r.farkas_row_upper[i], C.md.row_up[i], -1)) return false;
    double residual = 0;
    for (size_t j = 0; j < n; ++j) {
        double z = C.zcol(j, r.farkas_row_lower, r.farkas_row_upper, r.farkas_col_lower, r.farkas_col_upper);
        if (!std::isfinite(z)) return false;
        residual = std::max(residual, std::abs(z));
        if (!term(r.farkas_col_lower[j], C.lo[j], 1) || !term(r.farkas_col_upper[j], C.up[j], -1)) return false;
        if (z) {
            double b = C.need_bound(j, z);
            if (!std::isfinite(b)) return false;
            double t = z * b; if (!std::isfinite(t)) return false;
            terms.push_back(t);
        }
    }
    std::vector<double> ab; for (double t : terms) ab.push_back(std::abs(t));
    double margin = exact_sum(terms) / (1 + exact_sum(ab));
    return std::isfinite(norm) && std::isfinite(margin) && std::isfinite(residual) &&
           std::abs(norm - 1) <= 1e-8 && residual <= 1e-12 && margin > 1e-8;
}

static void repair_farkas(StrictCtx& C, Result& r, int max_iter = 400) {
    const size_t n = C.n, m = C.m; const double INF = std::numeric_limits<double>::infinity();
    auto& rl = r.farkas_row_lower; auto& ru = r.farkas_row_upper;
    auto& cl = r.farkas_col_lower; auto& cu = r.farkas_col_upper;
    std::vector<std::vector<int>> rowcols(m);
    for (size_t j = 0; j < n; ++j) for (const Entry& e : C.md.cols[j]) rowcols[e.index].push_back((int)j);
    auto bad = [&](size_t j) {
        double z = C.zcol(j, rl, ru, cl, cu);
        if (std::isnan(z)) return 1e300;
        if (!z) return 0.0;
        return std::isfinite(C.need_bound(j, z)) ? 0.0 : std::abs(z);
    };
    auto score = [&](const std::vector<int>& idx) {
        int cnt = 0; double sum = 0;
        for (int j : idx) { double b = bad(j); if (b) { ++cnt; sum += b; } }
        return std::make_pair(cnt, sum);
    };
    static const double ks[] = {1, 2, 3, 0.5, 1.5, 4, 8, 0.25, 1.01, 0.99, 2.5, 6};
    const long long kWorkCap = 20000000;  // deterministic budget (column-term evaluations), not wall clock
    for (int it = 0; it < max_iter; ++it) {
        bool prog = false, any = false;
        for (size_t j = 0; j < n; ++j) {
            if (C.work > kWorkCap) return;
            if (!bad(j)) continue;
            any = true;
            double z = C.zcol(j, rl, ru, cl, cu);
            if (!z || !bad(j)) continue;
            const double l = C.lo[j], u = C.up[j];
            bool done = false;
            for (double k : ks) {
                if (z > 0 && u < INF && cu[j] >= k * z) {
                    double old = cu[j]; cu[j] -= k * z;
                    if (bad(j) == 0) { prog = true; break; }
                    cu[j] = old;
                }
                if (z < 0 && l > -INF && cl[j] >= k * -z) {
                    double old = cl[j]; cl[j] = old - k * -z;
                    if (bad(j) == 0) { prog = true; break; }
                    cl[j] = old;
                }
                for (const Entry& e : C.md.cols[j]) {
                    int ri = e.index; double a = e.value;
                    if (a == 0) continue;
                    double delta = -k * z / a, rlo = C.md.row_lo[ri], rup = C.md.row_up[ri];
                    if (delta > 0 && !(ru[ri] > 0 || rup < INF)) continue;
                    if (delta < 0 && !(rl[ri] > 0 || rlo > -INF)) continue;
                    if (delta > 0 && rl[ri] > 0 && rl[ri] < delta) continue;
                    if (delta < 0 && ru[ri] > 0 && ru[ri] < -delta) continue;
                    auto before = score(rowcols[ri]);
                    double sl = rl[ri], su = ru[ri], ycur = su - sl;
                    if (delta > 0) { if (rl[ri] > 0) rl[ri] -= delta; else ru[ri] += delta; }
                    else { if (ru[ri] > 0) ru[ri] += delta; else rl[ri] -= delta; }
                    auto after = score(rowcols[ri]);
                    if ((ru[ri] - rl[ri]) == ycur || !(after < before)) { rl[ri] = sl; ru[ri] = su; }
                    else { prog = true; done = true; break; }
                }
                if (done) break;
            }
        }
        if (!any || !prog) break;
    }
}



// Bounds on the repair LP (the auxiliary solve_lp call only): its time limit is min(20 s per call, 60 s cumulative
// per thread, time left to the caller's deadline), and a thread-local guard stops the auxiliary solve from starting
// another repair LP. NOT deadline-bounded: StrictCtx construction, the local repair_farkas pass (bounded by its own
// work counter, not by the clock), auxiliary model construction, and the 60000-column skip test, which all run
// before or outside that budget. The budget is wall-clock based, so under extreme load a repair that would have
// fit may be skipped; the original certificate is then kept and flagged unresolved.
static thread_local int g_resolve_depth = 0;
static thread_local double g_resolve_spent_s = 0.0;
static constexpr double kResolvePerCallS = 20.0, kResolveTotalS = 60.0;
static constexpr size_t kResolveMaxCols = 60000;
static double resolve_budget(std::chrono::steady_clock::time_point deadline) {
    if (g_resolve_depth > 0) return 0.0;
    double b = std::min(kResolvePerCallS, kResolveTotalS - g_resolve_spent_s);
    if (deadline != std::chrono::steady_clock::time_point::max()) {
        double rem = std::chrono::duration<double>(deadline - std::chrono::steady_clock::now()).count();
        b = std::min(b, rem);
    }
    return b > 0.05 ? b : 0.0;
}
struct ResolveScope {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    ResolveScope() { ++g_resolve_depth; }
    ~ResolveScope() { --g_resolve_depth; g_resolve_spent_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }
};
// Proximity L1 re-solve of the Farkas multipliers on the certificate's own LP (engine simplex only).
// Finds v >= 0 nearest (L1) to the normalised engine multipliers such that the stationarity of every column
// lies in a small window whose sign matches a finite needed bound, free columns are exactly 0, sum v = 1 and
// the bound contradiction keeps at least half of the original. Strict acceptance is decided afterwards by the
// strict port and the internal verifier, never by this routine.
static bool lp_resolve_farkas(StrictCtx& C, Result& r, bool relax, double time_s) {
    const size_t n = C.n, m = C.m; const double INF = std::numeric_limits<double>::infinity();
    const size_t N = 2 * m + 2 * n;
    if (N == 0 || 3 * N + n > kResolveMaxCols || !(time_s > 0)) return false;
    const double S = 1e6, ZLO = 3e-8, ZHI = 3e-7;
    std::vector<double> v0(N);
    for (size_t i = 0; i < m; ++i) { v0[i] = r.farkas_row_lower[i]; v0[m + i] = r.farkas_row_upper[i]; }
    for (size_t j = 0; j < n; ++j) { v0[2 * m + j] = r.farkas_col_lower[j]; v0[2 * m + n + j] = r.farkas_col_upper[j]; }
    std::vector<double> tmp = v0; double tot = exact_sum(tmp);
    if (!std::isfinite(tot) || !(tot > 0)) return false;
    for (double& x : v0) { if (!std::isfinite(x)) return false; x = std::max(x, 0.0) / tot * S; }
    // bound contradiction coefficient per variable
    std::vector<double> tv(N, 0.0);
    auto fin = [&](double x) { return std::isfinite(x) ? x : 0.0; };
    for (size_t i = 0; i < m; ++i) { tv[i] = fin(C.md.row_lo[i]); tv[m + i] = -fin(C.md.row_up[i]); }
    for (size_t j = 0; j < n; ++j) { tv[2 * m + j] = fin(C.lo[j]); tv[2 * m + n + j] = -fin(C.up[j]); }
    double t0 = 0; for (size_t k = 0; k < N; ++k) t0 += tv[k] * v0[k];
    // row entries of the checker-aggregated matrix
    std::vector<std::vector<Entry>> rowent(m);
    for (size_t j = 0; j < n; ++j) for (const Entry& e : C.pyagg[j]) if (e.value != 0) rowent[e.index].push_back({(int)j, e.value});
    const size_t R = N + 1 + 2 * n + 1;
    Model M; M.name = "farkas_resolve";
    M.row_names.assign(R, std::string()); M.col_names.assign(3 * N, std::string());
    M.cols.assign(3 * N + n, {}); M.cost.assign(3 * N + n, 0.0); M.col_lo.assign(3 * N + n, 0.0); M.col_up.assign(3 * N + n, INF);
    M.col_names.assign(3 * N + n, std::string());
    for (size_t k = N; k < 3 * N; ++k) M.cost[k] = 1.0;
    for (size_t j = 0; j < n; ++j) { M.cost[3 * N + j] = 1e3; M.col_up[3 * N + j] = ZLO; }  // elastic lower edge
    M.row_lo.assign(R, 0.0); M.row_up.assign(R, 0.0);
    for (size_t k = 0; k < N; ++k) { M.row_lo[k] = M.row_up[k] = v0[k]; }
    M.row_lo[N] = M.row_up[N] = S;
    for (size_t j = 0; j < n; ++j) {
        const bool fl = std::isfinite(C.lo[j]), fu = std::isfinite(C.up[j]);
                bool nocol = true; for (const Entry& e : C.pyagg[j]) if (e.value != 0) { nocol = false; break; }
        double zl, zu;
        if (fl && fu) { zl = -ZHI; zu = ZHI; }
        else if (fl) { zl = nocol ? 0.0 : ZLO; zu = ZHI; }
        else if (fu) { zl = -ZHI; zu = nocol ? 0.0 : -ZLO; }
        else { zl = 0; zu = 0; }
        M.row_lo[N + 1 + j] = zl; M.row_up[N + 1 + j] = INF;  // z + e >= zl
        M.row_lo[N + 1 + n + j] = -INF; M.row_up[N + 1 + n + j] = zu;  // z <= zu
    }
    M.row_lo[R - 1] = 0.5 * t0; M.row_up[R - 1] = INF;
    for (size_t k = 0; k < N; ++k) {
        auto& col = M.cols[k];
        col.push_back({(int)k, 1.0}); col.push_back({(int)N, 1.0});
        std::vector<Entry> zc;  // (z column, coefficient), ascending
        if (k < m) { for (const Entry& e : rowent[k]) zc.push_back({e.index, -e.value}); }
        else if (k < 2 * m) { for (const Entry& e : rowent[k - m]) zc.push_back({e.index, e.value}); }
        else if (k < 2 * m + n) zc.push_back({(int)(k - 2 * m), -1.0});
        else zc.push_back({(int)(k - 2 * m - n), 1.0});
        for (const Entry& e : zc) col.push_back({(int)(N + 1 + e.index), e.value});
        for (const Entry& e : zc) col.push_back({(int)(N + 1 + n + e.index), e.value});
        if (tv[k] != 0) col.push_back({(int)(R - 1), tv[k]});
    }
    for (size_t j = 0; j < n; ++j) M.cols[3 * N + j].push_back({(int)(N + 1 + j), 1.0});
    for (size_t k = 0; k < N; ++k) { M.cols[N + k].push_back({(int)k, -1.0}); M.cols[2 * N + k].push_back({(int)k, 1.0}); }
    // v bounds: a multiplier on an infinite side is 0
    for (size_t i = 0; i < m; ++i) {
        if (!std::isfinite(C.md.row_lo[i])) M.col_up[i] = 0;
        if (!std::isfinite(C.md.row_up[i])) M.col_up[m + i] = 0;
    }
    for (size_t j = 0; j < n; ++j) {
        if (!std::isfinite(C.lo[j])) M.col_up[2 * m + j] = 0;
        if (!std::isfinite(C.up[j])) M.col_up[2 * m + n + j] = 0;
    }
    Result q;
    { ResolveScope guard; q = solve_lp(M, M.col_lo, M.col_up, nullptr, time_s); }
    if (q.status != Status::Optimal || q.x.size() < N) return false;
    std::vector<double> v(N); double sum = 0; std::vector<double> sv;
    for (size_t k = 0; k < N; ++k) { double x = q.x[k] / S; if (!(x > 0) || !std::isfinite(x)) x = 0; v[k] = x; sv.push_back(x); }
    sum = exact_sum(sv); if (!(sum > 0) || !std::isfinite(sum)) return false;
    for (double& x : v) x /= sum;
    for (size_t i = 0; i < m; ++i) { r.farkas_row_lower[i] = v[i]; r.farkas_row_upper[i] = v[m + i]; }
    for (size_t j = 0; j < n; ++j) { r.farkas_col_lower[j] = v[2 * m + j]; r.farkas_col_upper[j] = v[2 * m + n + j]; }
    return true;
}

bool verify_nonoptimal_certificate(const Model& md, const std::vector<double>& lo,
                                  const std::vector<double>& up, Result& r,
                                  std::chrono::steady_clock::time_point deadline) {
    // Keep the original certificate if the strict contract already holds or cannot be reached.
    const bool candidate = r.status == Status::Infeasible && r.farkas_row_lower.size() == md.row_names.size() &&
                           r.farkas_col_lower.size() == md.cols.size() && r.farkas_row_upper.size() == md.row_names.size() &&
                           r.farkas_col_upper.size() == md.cols.size() && lo.size() == md.cols.size() && up.size() == md.cols.size();
    if (!candidate) return verify_impl(md, lo, up, r);
    Result orig = r;
    bool ok = verify_impl(md, lo, up, r);
    StrictCtx C(md, lo, up);
    if (!ok) {
        // The solver-internal verifier rejected the multipliers (typically roundoff stationarity on
        // an unbounded column). Try the repair; adopt ONLY if the strict port AND the internal
        // verifier both accept the repaired multipliers. Otherwise return the original failure.
        Result rep = orig;
        repair_farkas(C, rep);
        if (strict_farkas_ok(C, rep)) {
            Result chk = rep;
            if (verify_impl(md, lo, up, chk)) { r = chk; return true; }
        }
        { Result lp = orig;
          if (lp_resolve_farkas(C, lp, false, resolve_budget(deadline)) && strict_farkas_ok(C, lp)) {
              Result chk = lp; if (verify_impl(md, lo, up, chk)) { r = chk; return true; } } }
        return false;
    }
    if (strict_farkas_ok(C, orig)) return true;
    Result rep = orig;
    repair_farkas(C, rep);
    if (strict_farkas_ok(C, rep)) {
        Result chk = rep;
        if (verify_impl(md, lo, up, chk)) { r = chk; return true; }
    }
    { Result lp = orig;
      if (lp_resolve_farkas(C, lp, false, resolve_budget(deadline)) && strict_farkas_ok(C, lp)) {
          Result chk = lp; if (verify_impl(md, lo, up, chk)) { r = chk; r.message += " [strict_farkas_repair=lp_resolve]"; return true; } } }
    // Original certificate kept: it passed the solver-internal verifier but NOT the strict contract.
    // Marked explicitly; certificate_verified keeps its old (solver-internal) meaning.
    r.message += " [strict_farkas_contract=unresolved]";
    return true;
}
