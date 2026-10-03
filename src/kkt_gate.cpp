// KKT-gated optimal with one equilibrated retry (feature/kkt-gated-optimal).
// An LP result is reported Optimal only when an independent check on the ORIGINAL model agrees:
// primal feasibility, dual feasibility (reduced costs judged against the cost scale of their own
// column) and complementarity at the returned basis. Duals come from a fresh LU of the returned basis,
// not from solver state. If the solve ends numerical_failure or the check fails, one retry runs on a
// Ruiz-equilibrated copy whose scale factors are all powers of two (exact in binary FP), mapped back
// exactly, and re-checked on the original model. Otherwise the status is numerical_failure.
// Infeasible / unbounded / time-limit results pass through untouched. C++17, standard library only.
#include <algorithm>
#include <chrono>
#include <cmath>

#include "taral.hpp"

namespace {
constexpr double kPrimalGate = 1e-7;  // same bar as Simplex::finish
constexpr double kDualGate = 1e-7;    // relative to the column's own cost scale (>= 1)

double pow2_round(double f) {
    if (!(f > 0) || !std::isfinite(f)) return 1.0;
    return std::ldexp(1.0, int(std::lround(std::log2(f))));
}
}  // namespace

KktReport kkt_check(const Model& md, const std::vector<double>& clo, const std::vector<double>& cup, const Result& r) {
    KktReport k;
    const int n = int(md.col_names.size()), m = int(md.row_names.size());
    auto fail = [&](const std::string& why) { k.ok = false; k.msg = why; return k; };
    if (int(r.x.size()) != n) return fail("no primal point");
    if (int(r.basis.size()) != n + m) return fail("no basis to certify");
    for (double v : r.x) if (!std::isfinite(v)) return fail("nonfinite primal point");
    const double sgn = md.maximize ? -1.0 : 1.0;  // work in minimisation form
    std::vector<double> act(m, 0.0);
    for (int j = 0; j < n; ++j) for (const Entry& e : md.cols[j]) act[e.index] += e.value * r.x[j];
    // primal feasibility on the original model
    for (int i = 0; i < m; ++i) {
        double s = 1 + std::abs(std::isfinite(md.row_lo[i]) ? md.row_lo[i] : md.row_up[i]);
        k.primal = std::max({k.primal, (md.row_lo[i] - act[i]) / s, (act[i] - md.row_up[i]) / s});
    }
    for (int j = 0; j < n; ++j) k.primal = std::max({k.primal, clo[j] - r.x[j], r.x[j] - cup[j]});
    if (k.primal > kPrimalGate) return fail("primal infeasible by " + std::to_string(k.primal));
    // duals from a fresh factorisation of the returned basis
    std::vector<int> head;
    for (int j = 0; j < n + m; ++j) if (r.basis[j] == 0) head.push_back(j);
    if (int(head.size()) != m) return fail("basis size " + std::to_string(head.size()) + " != rows " + std::to_string(m));
    std::vector<std::vector<Entry>> cols(m);
    std::vector<double> cb(m, 0.0);
    for (int p = 0; p < m; ++p) {
        int j = head[p];
        if (j < n) { cols[p] = md.cols[j]; cb[p] = sgn * md.cost[j]; }
        else cols[p].push_back({j - n, -1.0});
    }
    SparseLU lu;
    if (!lu.factor(m, cols)) return fail("returned basis is singular");
    std::vector<double> y(m, 0.0);
    lu.btran(cb, y);
    for (double v : y) if (!std::isfinite(v)) return fail("nonfinite duals");
    auto value = [&](int j) { return j < n ? r.x[j] : act[j - n]; };
    auto lower = [&](int j) { return j < n ? clo[j] : md.row_lo[j - n]; };
    auto upper = [&](int j) { return j < n ? cup[j] : md.row_up[j - n]; };
    for (int j = 0; j < n + m; ++j) {
        double d, scale;
        if (j < n) {
            double rt = 0, ay = 0;
            for (const Entry& e : md.cols[j]) { ay += e.value * y[e.index]; rt += std::abs(e.value * y[e.index]); }
            d = sgn * md.cost[j] - ay;
            scale = std::max({1.0, std::abs(md.cost[j]), rt});
        } else {
            d = y[j - n];
            scale = std::max(1.0, std::abs(y[j - n]));
        }
        double tol = kDualGate * scale;
        double lo = lower(j), up = upper(j), x = value(j);
        double viol;
        if (r.basis[j] == 0) viol = std::abs(d);  // basic: reduced cost must vanish
        else {
            bool atl = std::isfinite(lo) && std::abs(x - lo) <= kPrimalGate * (1 + std::abs(lo));
            bool atu = std::isfinite(up) && std::abs(x - up) <= kPrimalGate * (1 + std::abs(up));
            if (atl && atu) viol = 0;                   // fixed: any sign
            else if (atl) viol = std::max(0.0, -d);     // at lower: d >= 0
            else if (atu) viol = std::max(0.0, d);      // at upper: d <= 0
            else if (!std::isfinite(lo) && !std::isfinite(up) && std::abs(x) <= kPrimalGate) viol = std::abs(d);  // free at 0
            else return fail("nonbasic variable " + std::to_string(j) + " is off its bound (complementarity)");
        }
        k.dual = std::max(k.dual, viol / scale);
        if (viol > tol) return fail("dual infeasible by " + std::to_string(viol) + " (scale " + std::to_string(scale) + ") at variable " + std::to_string(j));
    }
    k.ok = true;
    return k;
}

namespace {
// The primal simplex gets kPrimalShare of the time; if it returns time_limit, the dual simplex (a different path,
// which solves DFL001 where the primal stalls) gets everything that is left.
constexpr double kPrimalShare = 0.2;

// Work-based routing: the primal simplex may use kPrimalShare of the cumulative work budget.
// `limit` is the remaining wall-clock budget. Without --time-limit the caller passes 1e9 (effectively none), so the
// work cap alone decides the outcome; with an explicit --time-limit the wall clock is enforced by the engines and a
// time_limit status means the outcome depended on the clock (documented in docs/deterministic_solve.md).
Result run_one_work(const Model& md, const std::vector<double>& lo, const std::vector<double>& up,
                    const std::vector<char>* warm, double limit, bool use_dual, bool fallback) {
    WorkBudget& w = work_budget();
    const auto t0 = std::chrono::steady_clock::now();
    if (use_dual) return solve_lp_dual(md, lo, up, warm, limit);
    const long total_cap = w.cap;
    if (fallback) w.cap = w.used + std::max(1L, long(kPrimalShare * double(total_cap - w.used)));  // cap 0 means off, so at least 1
    Result r = solve_lp(md, lo, up, warm, limit);
    w.cap = total_cap;
    if (!fallback || r.status != Status::IterationLimit || r.message != "work limit" || w.used >= total_cap) return r;
    const double spent = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Result d = solve_lp_dual(md, lo, up, warm, limit - spent);
    d.iterations += r.iterations;
    d.message = "route=primal-work-budget-then-dual primal_iters=" + std::to_string(r.iterations) + "; dual simplex fallback: " + d.message;
    return d;
}

Result run_one(const Model& md, const std::vector<double>& lo, const std::vector<double>& up,
               const std::vector<char>* warm, double limit, bool use_dual, bool fallback) {
    if (work_budget().cap) return run_one_work(md, lo, up, warm, limit, use_dual, fallback);
    if (use_dual) return solve_lp_dual(md, lo, up, warm, limit);
    const auto t0 = std::chrono::steady_clock::now();
    Result r = solve_lp(md, lo, up, warm, fallback ? kPrimalShare * limit : limit);
    if (!fallback || r.status != Status::TimeLimit) return r;
    const double spent = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::fprintf(stderr, "primal simplex not finished after %.1f s, switching to dual simplex\n", spent);
    Result d = solve_lp_dual(md, lo, up, warm, limit - spent);
    d.iterations += r.iterations;
    d.message = "route=primal-budget-then-dual budget_s=" + std::to_string(kPrimalShare * limit) + " primal_iters=" +
                std::to_string(r.iterations) + " primal_stalled_s=" + std::to_string(spent) + "; dual simplex fallback: " + d.message;
    return d;
}
}  // namespace

WorkBudget& work_budget() {
    static WorkBudget w;
    return w;
}

Result solve_lp_gated(const Model& md, const std::vector<double>& lo, const std::vector<double>& up,
                      const std::vector<char>* warm, double time_limit_s, bool use_dual, bool primal_fallback) {
    auto t0 = std::chrono::steady_clock::now();
    auto left = [&] { return time_limit_s - std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    Result r = run_one(md, lo, up, warm, time_limit_s, use_dual, primal_fallback);
    if (r.status != Status::Optimal && r.status != Status::NumericalFailure) return r;
    std::string why;
    if (r.status == Status::Optimal) {
        KktReport k = kkt_check(md, lo, up, r);
        if (k.ok) return r;
        why = "KKT gate: " + k.msg;
    } else {
        why = "first solve: " + r.message;
    }
    if ((!work_budget().cap || time_limit_s < 1e8) && left() < 0.5) {  // an unlimited wall budget (work mode) never stops here
        r.status = Status::NumericalFailure;
        r.message = why + "; no time left for the equilibrated retry";
        return r;
    }
    // Ruiz equilibration of A, every factor rounded to a power of two.
    const int n = int(md.col_names.size()), m = int(md.row_names.size());
    std::vector<double> rs(m, 1.0), cs(n, 1.0);
    for (int it = 0; it < 12; ++it) {
        std::vector<double> rmax(m, 0.0), cmax(n, 0.0);
        for (int j = 0; j < n; ++j)
            for (const Entry& e : md.cols[j]) {
                double a = std::abs(e.value) * rs[e.index] * cs[j];
                rmax[e.index] = std::max(rmax[e.index], a);
                cmax[j] = std::max(cmax[j], a);
            }
        for (int i = 0; i < m; ++i) if (rmax[i] > 0) rs[i] /= std::sqrt(rmax[i]);
        for (int j = 0; j < n; ++j) if (cmax[j] > 0) cs[j] /= std::sqrt(cmax[j]);
    }
    for (double& v : rs) v = pow2_round(v);
    for (double& v : cs) v = pow2_round(v);
    Model sm = md;
    for (int j = 0; j < n; ++j) {
        for (Entry& e : sm.cols[j]) e.value = e.value * rs[e.index] * cs[j];  // exact: powers of two
        sm.cost[j] = md.cost[j] * cs[j];
    }
    for (int i = 0; i < m; ++i) { sm.row_lo[i] = md.row_lo[i] * rs[i]; sm.row_up[i] = md.row_up[i] * rs[i]; }
    std::vector<double> slo(n), sup(n);
    for (int j = 0; j < n; ++j) { slo[j] = lo[j] / cs[j]; sup[j] = up[j] / cs[j]; sm.col_lo[j] = slo[j]; sm.col_up[j] = sup[j]; }
    Result s = run_one(sm, slo, sup, nullptr, left(), use_dual, false);
    if (s.status != Status::Optimal) {
        // a retry that proves infeasible/unbounded is not trusted over the first answer: report the failure
        Result f;
        // a retry stopped by a limit keeps that limit status (work cap or clock); anything else is a failure
        f.status = (s.status == Status::TimeLimit) ? Status::TimeLimit
                   : (s.status == Status::IterationLimit && work_budget().cap) ? Status::IterationLimit : Status::NumericalFailure;
        f.iterations = r.iterations + s.iterations;
        f.message = why + "; equilibrated retry ended " + status_name(s.status) + (s.message.empty() ? "" : " (" + s.message + ")");
        return f;
    }
    Result o = s;  // map back exactly
    o.x.assign(n, 0.0);
    for (int j = 0; j < n; ++j) o.x[j] = s.x[j] * cs[j];
    o.objective = md.obj_const;
    for (int j = 0; j < n; ++j) o.objective += md.cost[j] * o.x[j];
    o.iterations = r.iterations + s.iterations;
    KktReport k = kkt_check(md, lo, up, o);
    if (!k.ok) {
        o.status = Status::NumericalFailure;
        o.message = why + "; equilibrated retry also failed the KKT gate: " + k.msg;
        return o;
    }
    o.message = (o.message.empty() ? "" : o.message + "; ") + "recovered by equilibrated retry (" + why + ")";
    return o;
}
