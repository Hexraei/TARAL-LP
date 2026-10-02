// Original-coordinate LP proof checker. No presolve data or basis flags are trusted.
#include <algorithm>
#include <cmath>
#include "taral.hpp"

bool verify_nonoptimal_certificate(const Model& md, const std::vector<double>& lo,
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
