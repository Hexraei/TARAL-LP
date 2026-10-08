// PROTOTYPE (never ships in src/ as-is): convex MIQP by branch and bound over QP relaxations.
// Each node relaxation is solved by ipm_solve on the model with the node's column bounds.
// Scope and honesty limits:
//  - minimization only (maximize + quadratic objective returns "unsupported");
//  - convex diagonal Q only: every qobj entry must have row == col and value >= -1e-12
//    (off-diagonal or negative diagonal returns "unsupported" - general convexity witnessing
//    is a separate backlog item, not silently approximated here);
//  - no presolve, propagation, cuts or warm starts; small cases only.
// Tolerances mirror src/milp.cpp: kIntTol = kFeasTol = kGapTol = 1e-6 (same conventions, same provenance).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <queue>
#include <string>
#include <vector>
#include "taral.hpp"
#include "ipm.hpp"
#include "miqp.hpp"

namespace miqp_proto {
constexpr double kIntTol = 1e-6;
constexpr double kGapTol = 1e-6;
constexpr double kInf = std::numeric_limits<double>::infinity();



static double quad_objective_of(const Model& md, const std::vector<double>& x) {
    double z = md.obj_const;
    for (size_t j = 0; j < x.size(); ++j) z += md.cost[j] * x[j];
    for (const QEntry& q : md.qobj) {
        double t = q.value * x[q.row] * x[q.col];
        z += (q.row == q.col) ? 0.5 * t : t;  // 0.5 x'Qx with lower triangle stored once
    }
    return z;
}

static double violation_lin(const Model& md, const std::vector<double>& x) {
    int n = int(md.cols.size()), m = int(md.row_lo.size());
    std::vector<double> act(m, 0.0);
    double worst = 0;
    auto check = [&](double v, double lo, double up) {
        if (std::isfinite(lo)) worst = std::max(worst, (lo - v) / (1 + std::abs(lo)));
        if (std::isfinite(up)) worst = std::max(worst, (v - up) / (1 + std::abs(up)));
    };
    for (int j = 0; j < n; ++j) {
        if (!std::isfinite(x[j])) return kInf;
        for (const Entry& e : md.cols[j]) act[e.index] += e.value * x[j];
        check(x[j], md.col_lo[j], md.col_up[j]);
        if (md.is_int[j]) worst = std::max(worst, std::abs(x[j] - std::round(x[j])));
    }
    for (int i = 0; i < m; ++i) check(act[i], md.row_lo[i], md.row_up[i]);
    return worst;
}

struct Node {
    double bound;
    long seq;
    std::vector<double> lo, up;
};
struct Worse { bool operator()(const Node& a, const Node& b) const { return a.bound > b.bound; } };

static std::string convex_diag_check(const Model& md) {
    if (md.maximize) return "prototype supports minimization only";
    for (const QEntry& q : md.qobj) {
        if (q.row != q.col) return "prototype supports diagonal Q only (off-diagonal term present)";
        if (q.value < -1e-12) return "prototype supports convex Q only (negative diagonal term)";
    }
    return "";
}

MiqpResult solve_miqp(const Model& md, double time_limit, long node_limit) {
    MiqpResult res;
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    std::string why = convex_diag_check(md);
    if (!why.empty()) { res.status = "unsupported"; res.message = why; return res; }
    const int n = int(md.cols.size());
    double z = kInf;
    std::priority_queue<Node, std::vector<Node>, Worse> open;
    long seq = 0;
    open.push(Node{-kInf, seq++, md.col_lo, md.col_up});
    double bound = -kInf;
    bool hit_limit = false;
    while (!open.empty()) {
        Node cur = open.top(); open.pop();
        if (cur.bound > z - kGapTol * std::max(1.0, std::abs(z))) continue;  // pruned by incumbent
        if (elapsed() > time_limit) { hit_limit = true; open.push(std::move(cur)); break; }
        if (node_limit >= 0 && res.nodes >= node_limit) { hit_limit = true; open.push(std::move(cur)); break; }
        Model nm = md;
        nm.col_lo = cur.lo; nm.col_up = cur.up;
        nm.is_int.assign(nm.is_int.size(), 0);  // relaxation: continuous columns (ipm refuses integer columns)
        IpmOptions opt; opt.time_limit = std::max(0.001, time_limit - elapsed());
        IpmResult r = ipm_solve(nm, opt);
        ++res.nodes;
        if (r.status == IpmStatus::Infeasible) continue;  // prune
        if (r.status != IpmStatus::Optimal) { ++res.unresolved; continue; }  // unbounded/dual-infeasible/failure: no valid bound, do not prune or branch
        double nb = r.objective;  // QP relaxation optimum = node lower bound (min sense)
        bound = open.empty() ? nb : std::min(nb, open.top().bound);
         if (nb > z - kGapTol * std::max(1.0, std::abs(z))) continue;  // bound-pruned
        // incumbent offer
        if (r.x.size() == static_cast<size_t>(n)) {
            std::vector<double> xr = r.x;
            for (int j = 0; j < n; ++j) if (md.is_int[j]) xr[j] = std::round(xr[j]);
            for (const std::vector<double>* c : {&xr, &r.x}) {
                if (violation_lin(md, *c) > kIntTol) continue;
                double v = quad_objective_of(md, *c);
                if (v < z) { z = v; res.x = *c; res.has_solution = true; }
                break;
            }
        }
        // branch on the most fractional integer column
        int bj = -1; double bf = 0;
        for (int j = 0; j < n; ++j) {
            if (!md.is_int[j]) continue;
            double f = std::abs(r.x[j] - std::round(r.x[j]));
            if (f > kIntTol && f > bf) { bf = f; bj = j; }
        }
        if (bj < 0) {
            // near-integer relaxation optimum: rounding alone can move the objective by more
            // than the gap tolerance (observed: case409, 1.3e-6 suboptimal incumbent at
            // f <= kIntTol). Polish: fix integer columns at their rounded values, resolve the
            // restricted continuous QP, and offer that exact point. Fall back to the rounded
            // point if the restricted solve fails.
            Model pm = nm;
            std::vector<double> xr = r.x;
            for (int j = 0; j < n; ++j) if (md.is_int[j]) {
                xr[j] = std::min(pm.col_up[j], std::max(pm.col_lo[j], std::round(r.x[j])));
                pm.col_lo[j] = pm.col_up[j] = xr[j];
            }
            IpmOptions popt; popt.time_limit = std::max(0.001, time_limit - elapsed());
            IpmResult pr = ipm_solve(pm, popt);
            ++res.nodes;
            if (pr.status == IpmStatus::Optimal && pr.x.size() == static_cast<size_t>(n)
                && violation_lin(md, pr.x) <= kIntTol) {
                double v = quad_objective_of(md, pr.x);
                if (v < z) { z = v; res.x = pr.x; res.has_solution = true; }
            } else if (violation_lin(md, xr) <= kIntTol) {
                double v = quad_objective_of(md, xr);
                if (v < z) { z = v; res.x = xr; res.has_solution = true; }
            }
            continue;
        }
        double v0 = r.x[bj];
        Node dn{nb, seq++, cur.lo, cur.up}; dn.up[bj] = std::floor(v0);
        Node up{nb, seq++, cur.lo, cur.up}; up.lo[bj] = std::ceil(v0);
        if (dn.up[bj] >= dn.lo[bj]) open.push(std::move(dn));
        if (up.lo[bj] <= up.up[bj]) open.push(std::move(up));
    }
    res.objective = z;
    double rem = open.empty() ? kInf : open.top().bound;
    res.best_bound = std::min(bound, rem);
    if (res.best_bound == -kInf && res.nodes > 0) res.best_bound = bound;
    if (hit_limit) { res.status = "time_limit"; res.message = "cap reached with open nodes"; }
    else if (!open.empty() && res.unresolved > 0) { res.status = "time_limit"; res.message = "unresolved nodes remain"; }
    else if (res.has_solution) res.status = "optimal";
    else if (res.unresolved == 0) res.status = "infeasible";
    else { res.status = "numerical_failure"; res.message = "all nodes unresolved"; }
    return res;
}
}  // namespace miqp_proto
