// Verified presolve: a small set of exact-structure reductions with a replayable log and an original-space audit.
// Log field meaning per op (PresolveOp: row, col, a, v1..v4):
//   round_int_bounds  col, v1/v2 = old lo/up, v3/v4 = new lo/up (ceil/floor with 1e-9 relative slack)
//   fix_col           col, v1 = value (lo == up in the current state); rows' bounds and the objective constant shift
//   empty_row         row, v1/v2 = current row bounds; 0 lies inside them (within 1e-9 relative)
//   singleton_row     row, col, a = coefficient, v1/v2 = current row bounds, v3/v4 = column bounds afterwards
//   redundant_row     row, v1/v2 = row bounds, v3/v4 = min/max activity over current column bounds (finite)
//   fix_empty_col     col, v1 = value; the column is in no remaining row and the cost sign picks the bound
// Infeasibility is only reported with the reason that the replayer re-derives (same tolerances).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "taral.hpp"

namespace {
constexpr double kRel = 1e-9;
double sc(double b) { return 1 + (std::isfinite(b) ? std::abs(b) : 0); }
// Strict inward integer rounding, no outward snap: ceil(lb), floor(ub). An interval that contains no integer
// ([4e-7,8e-7], [1.0000004,1.0000008]) is infeasible exactly; nothing is approximated.
double int_up(double l) { return std::ceil(l); }
double int_dn(double u) { return std::floor(u); }
// Audit scale for column bounds: relative up to 1e3, absolute (1e-6 * 1e3) beyond, so large bounds stay strict.
double scb(double b) { return 1 + (std::isfinite(b) ? std::min(std::abs(b), 1e3) : 0); }
}  // namespace

PresolveResult presolve_model(const Model& orig, std::chrono::steady_clock::time_point deadline, long test_timeout_after_ops) {
    PresolveResult P;
    long ticks = 0;
    auto expired = [&]() {
        if (P.timed_out) return true;
        if (test_timeout_after_ops >= 0 && (long)P.log.size() >= test_timeout_after_ops) return P.timed_out = true;
        if (std::chrono::steady_clock::now() > deadline) return P.timed_out = true;
        return false;
    };
    auto tick = [&]() { return (test_timeout_after_ops >= 0 || (++ticks & 255) == 0) && expired(); };
    const int n = (int)orig.cols.size(), m = (int)orig.row_lo.size();
    std::vector<double> lo = orig.col_lo, up = orig.col_up, rlo = orig.row_lo, rup = orig.row_up;
    double objc = orig.obj_const;
    std::vector<std::vector<Entry>> rows(m);  // (col, value)
    for (int j = 0; j < n; ++j)
        for (const Entry& e : orig.cols[j])
            if (e.value != 0) rows[e.index].push_back({j, e.value});
    std::vector<char> rrem(m, 0), crem(n, 0);
    P.fixed_value.assign(n, 0.0);
    auto fail = [&](const std::string& why) { P.infeasible = true; P.infeasible_reason = why; };

    if (expired()) { P.removed_col.assign(n, 0); P.fixed_value.assign(n, 0.0); return P; }
    for (int j = 0; j < n && !P.infeasible && !P.timed_out; ++j) {
        if (tick()) break;
        if (!orig.is_int.empty() && orig.is_int[j]) {
            double nl = std::isfinite(lo[j]) ? int_up(lo[j]) : lo[j];
            double nu = std::isfinite(up[j]) ? int_dn(up[j]) : up[j];
            if (nl != lo[j] || nu != up[j]) {
                P.log.push_back({"round_int_bounds", -1, j, 0, lo[j], up[j], nl, nu});
                lo[j] = nl, up[j] = nu;
            }
        }
        const bool jint = !orig.is_int.empty() && orig.is_int[j];
        if (jint ? lo[j] > up[j] : lo[j] > up[j] + kRel * sc(up[j])) fail("column " + std::to_string(j) + " has lower bound above upper bound");
    }
    bool changed = !P.infeasible && !P.timed_out;
    while (changed && P.passes < 50) {
        if (expired()) break;
        changed = false;
        ++P.passes;
        for (int j = 0; j < n && !P.infeasible && !P.timed_out; ++j) {  // fixed columns
            if (tick()) break;
            if (crem[j] || !(lo[j] == up[j]) || !std::isfinite(lo[j])) continue;
            double v = lo[j];
            P.log.push_back({"fix_col", -1, j, 0, v, 0, 0, 0});
            crem[j] = 1, P.fixed_value[j] = v, changed = true;
            for (const Entry& e : orig.cols[j]) {
                if (rrem[e.index]) continue;
                if (std::isfinite(rlo[e.index])) rlo[e.index] -= e.value * v;
                if (std::isfinite(rup[e.index])) rup[e.index] -= e.value * v;
            }
            objc += orig.cost[j] * v;
        }
        for (int i = 0; i < m && !P.infeasible && !P.timed_out; ++i) {
            if (tick()) break;
            if (rrem[i]) continue;
            int cnt = 0, last = -1;
            double lastv = 0;
            for (const Entry& e : rows[i])
                if (!crem[e.index]) ++cnt, last = e.index, lastv = e.value;
            if (cnt == 0) {
                if (rlo[i] > kRel * sc(rlo[i]) || rup[i] < -kRel * sc(rup[i])) { fail("row " + std::to_string(i) + " empty and violated"); break; }
                P.log.push_back({"empty_row", i, -1, 0, rlo[i], rup[i], 0, 0});
                rrem[i] = 1, changed = true;
                continue;
            }
            if (cnt == 1) {
                double l = rlo[i] / lastv, u = rup[i] / lastv;
                if (lastv < 0) std::swap(l, u);
                if (!orig.is_int.empty() && orig.is_int[last]) {
                    if (std::isfinite(l)) l = int_up(l);
                    if (std::isfinite(u)) u = int_dn(u);
                }
                double nl = std::max(lo[last], l), nu = std::min(up[last], u);
                const bool lint = !orig.is_int.empty() && orig.is_int[last];  // integer interval: exact comparison
                if (lint ? nl > nu : nl > nu + kRel * sc(nu)) { fail("singleton row " + std::to_string(i) + " conflicts with bounds of column " + std::to_string(last)); break; }
                if (nl > nu) nu = nl;  // equal within tolerance
                lo[last] = nl, up[last] = nu;
                P.log.push_back({"singleton_row", i, last, lastv, rlo[i], rup[i], nl, nu});
                rrem[i] = 1, changed = true;
                continue;
            }
            double mn = 0, mx = 0;
            for (const Entry& e : rows[i]) {
                if (crem[e.index]) continue;
                double a = e.value, l = lo[e.index], u = up[e.index];
                mn += a > 0 ? a * l : a * u;
                mx += a > 0 ? a * u : a * l;
            }
            if (std::isfinite(mn) && std::isfinite(mx)) {
                if (mn > rup[i] + kRel * sc(rup[i]) || mx < rlo[i] - kRel * sc(rlo[i])) {
                    fail("row " + std::to_string(i) + " activity range [" + std::to_string(mn) + "," + std::to_string(mx) + "] outside row bounds");
                    break;
                }
            }
            if (std::isfinite(mn) && std::isfinite(mx) && mn >= rlo[i] && mx <= rup[i]) {
                P.log.push_back({"redundant_row", i, -1, 0, rlo[i], rup[i], mn, mx});
                rrem[i] = 1, changed = true;
            } else if (!std::isfinite(rlo[i]) && std::isfinite(mx) && mx <= rup[i]) {
                P.log.push_back({"redundant_row", i, -1, 0, rlo[i], rup[i], mn, mx});
                rrem[i] = 1, changed = true;
            } else if (!std::isfinite(rup[i]) && std::isfinite(mn) && mn >= rlo[i]) {
                P.log.push_back({"redundant_row", i, -1, 0, rlo[i], rup[i], mn, mx});
                rrem[i] = 1, changed = true;
            }
        }
        for (int j = 0; j < n && !P.infeasible && !P.timed_out; ++j) {  // empty columns
            if (tick()) break;
            if (crem[j]) continue;
            bool any = false;
            for (const Entry& e : orig.cols[j])
                if (!rrem[e.index] && e.value != 0) { any = true; break; }
            if (any) continue;
            double c = orig.maximize ? -orig.cost[j] : orig.cost[j];
            double v;
            if (c > 0) v = lo[j]; else if (c < 0) v = up[j];
            else v = std::isfinite(lo[j]) ? lo[j] : (std::isfinite(up[j]) ? up[j] : 0.0);
            if (!std::isfinite(v)) continue;  // unbounded direction: leave to the solver
            P.log.push_back({"fix_empty_col", -1, j, 0, v, 0, 0, 0});
            lo[j] = up[j] = v;  // fix_col applies it next pass
            changed = true;
        }
    }
    // the loop only ends when no reduction applied (or the pass cap was hit); state is final
    P.removed_col = crem;
    if (P.infeasible || P.timed_out) return P;
    std::vector<int> rmap(m, -1), cmap(n, -1);
    Model& R = P.reduced;
    R.name = orig.name, R.maximize = orig.maximize, R.obj_const = objc;
    for (int i = 0; i < m; ++i)
        if (!rrem[i]) rmap[i] = (int)P.kept_rows.size(), P.kept_rows.push_back(i), R.row_names.push_back(orig.row_names[i]),
        R.row_lo.push_back(rlo[i]), R.row_up.push_back(rup[i]);
    for (int j = 0; j < n; ++j) {
        if (crem[j]) continue;
        cmap[j] = (int)P.kept_cols.size();
        P.kept_cols.push_back(j);
        R.col_names.push_back(orig.col_names[j]);
        R.cost.push_back(orig.cost[j]);
        R.col_lo.push_back(lo[j]);
        R.col_up.push_back(up[j]);
        R.is_int.push_back(orig.is_int.empty() ? 0 : orig.is_int[j]);
        std::vector<Entry> col;
        for (const Entry& e : orig.cols[j])
            if (rmap[e.index] >= 0) col.push_back({rmap[e.index], e.value});
        R.cols.push_back(col);
    }
    return P;
}

std::vector<double> presolve_expand(const Model& orig, const PresolveResult& p, const std::vector<double>& xr) {
    std::vector<double> x(orig.cols.size(), 0.0);
    for (size_t j = 0; j < x.size(); ++j) x[j] = p.removed_col[j] ? p.fixed_value[j] : 0.0;
    for (size_t k = 0; k < p.kept_cols.size(); ++k) x[p.kept_cols[k]] = xr[k];
    return x;
}

PresolveAudit presolve_audit(const Model& o, const std::vector<double>& x, double reported) {
    PresolveAudit a;
    const size_t n = o.cols.size(), m = o.row_lo.size();
    std::vector<double> ax(m, 0.0);
    double obj = o.obj_const;
    for (size_t j = 0; j < n; ++j) {
        obj += o.cost[j] * x[j];
        for (const Entry& e : o.cols[j]) ax[e.index] += e.value * x[j];
        a.max_bound_violation = std::max({a.max_bound_violation, (o.col_lo[j] - x[j]) / scb(o.col_lo[j]), (x[j] - o.col_up[j]) / scb(o.col_up[j])});
        if (!o.is_int.empty() && o.is_int[j]) a.max_int_violation = std::max(a.max_int_violation, std::abs(x[j] - std::round(x[j])));
    }
    for (size_t i = 0; i < m; ++i)
        a.max_row_violation = std::max({a.max_row_violation, (o.row_lo[i] - ax[i]) / sc(o.row_lo[i]), (ax[i] - o.row_up[i]) / sc(o.row_up[i])});
    a.objective = obj, a.reported_objective = reported, a.objective_diff = std::abs(obj - reported);
    a.ok = a.max_row_violation <= 1e-6 && a.max_bound_violation <= 1e-6 && a.max_int_violation <= 1e-6 &&
           a.objective_diff <= 1e-6 * (1 + std::abs(obj)) && std::isfinite(obj);
    return a;
}

static void jnum(std::FILE* f, double v) {
    if (std::isfinite(v)) std::fprintf(f, "%.17g", v);
    else std::fprintf(f, v > 0 ? "\"inf\"" : (v < 0 ? "\"-inf\"" : "null"));
}

bool write_presolve_log(const char* path, const Model& o, const PresolveResult& p, const PresolveAudit* au) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) return false;
    std::fprintf(f, "{\"original\": {\"rows\": %zu, \"cols\": %zu}, \"infeasible\": %s, \"timed_out\": %s, \"passes\": %d, \"tolerance_rel\": 1e-9",
                 o.row_lo.size(), o.cols.size(), p.infeasible ? "true" : "false", p.timed_out ? "true" : "false", p.passes);
    std::string r;
    for (char ch : p.infeasible_reason) r += (ch == '"' || ch == '\\') ? '_' : ch;
    std::fprintf(f, ", \"infeasible_reason\": \"%s\", \"ops\": [", r.c_str());
    for (size_t k = 0; k < p.log.size(); ++k) {
        const PresolveOp& op = p.log[k];
        std::fprintf(f, "%s{\"type\": \"%s\", \"row\": %d, \"col\": %d, \"a\": ", k ? ", " : "", op.type.c_str(), op.row, op.col);
        jnum(f, op.a);
        const double* vs[4] = {&op.v1, &op.v2, &op.v3, &op.v4};
        for (int q = 0; q < 4; ++q) { std::fprintf(f, ", \"v%d\": ", q + 1); jnum(f, *vs[q]); }
        std::fprintf(f, "}");
    }
    std::fprintf(f, "], \"kept_rows\": [");
    for (size_t k = 0; k < p.kept_rows.size(); ++k) std::fprintf(f, "%s%d", k ? ", " : "", p.kept_rows[k]);
    std::fprintf(f, "], \"kept_cols\": [");
    for (size_t k = 0; k < p.kept_cols.size(); ++k) std::fprintf(f, "%s%d", k ? ", " : "", p.kept_cols[k]);
    std::fprintf(f, "], \"reduced\": {\"rows\": %zu, \"cols\": %zu, \"obj_const\": ", p.kept_rows.size(), p.kept_cols.size());
    jnum(f, p.reduced.obj_const);
    std::fprintf(f, ", \"col_lo\": [");
    for (size_t k = 0; k < p.reduced.col_lo.size(); ++k) { if (k) std::fprintf(f, ", "); jnum(f, p.reduced.col_lo[k]); }
    std::fprintf(f, "], \"col_up\": [");
    for (size_t k = 0; k < p.reduced.col_up.size(); ++k) { if (k) std::fprintf(f, ", "); jnum(f, p.reduced.col_up[k]); }
    std::fprintf(f, "], \"row_lo\": [");
    for (size_t k = 0; k < p.reduced.row_lo.size(); ++k) { if (k) std::fprintf(f, ", "); jnum(f, p.reduced.row_lo[k]); }
    std::fprintf(f, "], \"row_up\": [");
    for (size_t k = 0; k < p.reduced.row_up.size(); ++k) { if (k) std::fprintf(f, ", "); jnum(f, p.reduced.row_up[k]); }
    std::fprintf(f, "]}");
    if (au) {
        std::fprintf(f, ", \"audit\": {\"ok\": %s, \"max_row_violation\": ", au->ok ? "true" : "false");
        jnum(f, au->max_row_violation);
        std::fprintf(f, ", \"max_bound_violation\": "); jnum(f, au->max_bound_violation);
        std::fprintf(f, ", \"max_int_violation\": "); jnum(f, au->max_int_violation);
        std::fprintf(f, ", \"objective\": "); jnum(f, au->objective);
        std::fprintf(f, ", \"reported_objective\": "); jnum(f, au->reported_objective);
        std::fprintf(f, ", \"objective_diff\": "); jnum(f, au->objective_diff);
        std::fprintf(f, "}");
    }
    std::fprintf(f, "}\n");
    bool ok = !std::ferror(f);
    if (std::fclose(f) != 0) ok = false;
    return ok;
}
