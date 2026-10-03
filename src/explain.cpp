// Verified LP infeasibility explanation: a row-irreducible subsystem plus a minimum-relaxation report.
// Scope and definitions are in docs/infeasibility_explanation.md. C++17 standard library only.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include "taral.hpp"

namespace {
using Clock = std::chrono::steady_clock;
constexpr double kWitnessTol = 1e-7;  // relative row/bound violation allowed for a "feasible" witness

// The model restricted to the rows in `keep` (column bounds untouched), all costs zero, row order preserved.
// map[k] = original row index of sub row k.
Model sub_model(const Model& m, const std::vector<char>& keep, std::vector<int>& map) {
    Model s;
    s.name = m.name;
    s.col_names = m.col_names;
    s.col_lo = m.col_lo;
    s.col_up = m.col_up;
    s.cost.assign(m.cols.size(), 0.0);
    std::vector<int> newidx(m.row_names.size(), -1);
    map.clear();
    for (size_t i = 0; i < m.row_names.size(); ++i)
        if (keep[i]) {
            newidx[i] = int(map.size());
            map.push_back(int(i));
            s.row_names.push_back(m.row_names[i]);
            s.row_lo.push_back(m.row_lo[i]);
            s.row_up.push_back(m.row_up[i]);
        }
    s.cols.resize(m.cols.size());
    for (size_t j = 0; j < m.cols.size(); ++j)
        for (const Entry& e : m.cols[j])
            if (newidx[e.index] >= 0) s.cols[j].push_back({newidx[e.index], e.value});
    s.is_int.assign(m.cols.size(), 0);
    return s;
}

// Worst relative violation of x against the rows in `keep` and all column bounds, long double sums.
double max_violation(const Model& m, const std::vector<char>& keep, const std::vector<double>& x) {
    double worst = 0;
    std::vector<long double> act(m.row_names.size(), 0);
    for (size_t j = 0; j < m.cols.size(); ++j) {
        if (!std::isfinite(x[j])) return kInf;
        if (std::isfinite(m.col_lo[j])) worst = std::max(worst, (m.col_lo[j] - x[j]) / (1 + std::abs(m.col_lo[j])));
        if (std::isfinite(m.col_up[j])) worst = std::max(worst, (x[j] - m.col_up[j]) / (1 + std::abs(m.col_up[j])));
        for (const Entry& e : m.cols[j]) act[e.index] += (long double)e.value * x[j];
    }
    for (size_t i = 0; i < m.row_names.size(); ++i) {
        if (!keep[i]) continue;
        if (std::isfinite(m.row_lo[i])) worst = std::max(worst, (double)((m.row_lo[i] - act[i]) / (1 + std::abs(m.row_lo[i]))));
        if (std::isfinite(m.row_up[i])) worst = std::max(worst, (double)((act[i] - m.row_up[i]) / (1 + std::abs(m.row_up[i]))));
    }
    return worst;
}

struct Outcome {
    enum Kind { Infeasible, Feasible, Unknown } kind = Unknown;
    Result r;                         // Infeasible: multipliers in original row indexing (full length)
    std::vector<double> x;            // Feasible: witness
    double violation = 0;
    std::string why;
};

Outcome test_rows(const Model& m, const std::vector<char>& keep, double budget_s) {
    Outcome o;
    std::vector<int> map;
    Model s = sub_model(m, keep, map);
    Result r = solve_lp_gated(s, s.col_lo, s.col_up, nullptr, budget_s, false, true);
    if (r.status == Status::Infeasible && r.certificate_verified &&
        r.farkas_row_lower.size() == map.size() && r.farkas_row_upper.size() == map.size()) {
        o.kind = Outcome::Infeasible;
        o.r = r;
        o.r.farkas_row_lower.assign(m.row_names.size(), 0.0);
        o.r.farkas_row_upper.assign(m.row_names.size(), 0.0);
        for (size_t k = 0; k < map.size(); ++k) {
            o.r.farkas_row_lower[map[k]] = r.farkas_row_lower[k];
            o.r.farkas_row_upper[map[k]] = r.farkas_row_upper[k];
        }
        return o;
    }
    if (r.status == Status::Optimal && r.x.size() == m.cols.size()) {
        double v = max_violation(m, keep, r.x);
        if (v <= kWitnessTol) {
            o.kind = Outcome::Feasible;
            o.x = r.x;
            o.violation = std::max(0.0, v);
            return o;
        }
        o.why = "solver point violates the subsystem by " + std::to_string(v);
        return o;
    }
    o.why = std::string("subsystem solve ended ") + status_name(r.status) + (r.message.empty() ? "" : ": " + r.message);
    return o;
}
}  // namespace

InfeasibilityExplanation explain_infeasibility(const Model& m, double time_limit_s) {
    const auto t0 = Clock::now();
    auto left = [&] { return time_limit_s - std::chrono::duration<double>(Clock::now() - t0).count(); };
    InfeasibilityExplanation e;
    const size_t nrows = m.row_names.size(), ncols = m.cols.size();
    long solves = 0;
    auto finish = [&](const char* status, const std::string& msg) {
        e.status = status;
        e.message = msg;
        e.lp_solves = solves;
        e.wall_s = std::chrono::duration<double>(Clock::now() - t0).count();
    };
    for (size_t j = 0; j < ncols; ++j)
        if (!(m.col_lo[j] <= m.col_up[j] + 0.0) && !std::isnan(m.col_lo[j]) && !std::isnan(m.col_up[j]))
            e.inconsistent_bound_columns.push_back(int(j));

    std::vector<char> all(nrows, 1);
    ++solves;
    Outcome root = test_rows(m, all, left());
    if (root.kind == Outcome::Feasible) { finish("relaxation_feasible", "the LP relaxation is feasible; nothing to explain"); return e; }
    if (root.kind == Outcome::Unknown) { finish("no_verified_proof", "no verified infeasibility proof for the full model: " + root.why); return e; }

    // Current certificate (multipliers in original row indexing) and its row support.
    Result cert = root.r;
    auto support = [&](const Result& r, const std::vector<char>& within) {
        std::vector<char> s(nrows, 0);
        for (size_t i = 0; i < nrows; ++i)
            if (within[i] && (r.farkas_row_lower[i] != 0 || r.farkas_row_upper[i] != 0)) s[i] = 1;
        return s;
    };
    std::vector<char> S = support(cert, all);
    // Deletion filter, least-weighted rows first.
    std::vector<int> order;
    for (size_t i = 0; i < nrows; ++i) if (S[i]) order.push_back(int(i));
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        return cert.farkas_row_lower[a] + cert.farkas_row_upper[a] < cert.farkas_row_lower[b] + cert.farkas_row_upper[b];
    });
    std::vector<std::vector<double>> witness(nrows);
    std::vector<double> wviol(nrows, 0.0);
    std::vector<char> unproven(nrows, 0);
    bool out_of_time = false;
    for (int r : order) {
        if (!S[r]) continue;
        if (left() <= 0) { unproven[r] = 1; out_of_time = true; continue; }
        std::vector<char> trial = S;
        trial[r] = 0;
        ++solves;
        Outcome o = test_rows(m, trial, left());
        if (o.kind == Outcome::Infeasible) {
            cert = o.r;
            S = support(cert, trial);  // rows without multiplier are not needed for this proof
        } else if (o.kind == Outcome::Feasible) {
            witness[r] = o.x;
            wviol[r] = o.violation;
        } else {
            unproven[r] = 1;
        }
    }
    // Final certificate must verify on the FULL model with the original column bounds.
    Result check = cert;
    check.status = Status::Infeasible;
    bool ok = verify_nonoptimal_certificate(m, m.col_lo, m.col_up, check);
    e.certificate_verified = ok;
    e.certificate_margin = check.certificate_margin;
    e.certificate_residual = check.certificate_residual;
    e.farkas_row_lower = cert.farkas_row_lower;
    e.farkas_row_upper = cert.farkas_row_upper;
    e.farkas_col_lower = cert.farkas_col_lower;
    e.farkas_col_upper = cert.farkas_col_upper;
    for (size_t i = 0; i < nrows; ++i)
        if (S[i]) {
            e.rows.push_back(int(i));
            e.witness.push_back(unproven[i] || witness[i].empty() ? std::vector<double>() : witness[i]);
            e.witness_violation.push_back(wviol[i]);
            if (unproven[i] || witness[i].empty()) e.unproven_rows.push_back(int(i));
        }
    for (size_t j = 0; j < ncols; ++j)
        if (j < e.farkas_col_lower.size() && (e.farkas_col_lower[j] != 0 || e.farkas_col_upper[j] != 0))
            e.bound_columns.push_back(int(j));
    if (!ok) {
        e.rows.clear(); e.witness.clear(); e.witness_violation.clear(); e.unproven_rows.clear(); e.bound_columns.clear();
        finish("no_verified_proof", "the reduced multiplier set failed the full-model proof check");
        return e;
    }

    // Minimum-relaxation report: weighted L1 side relaxation over ALL rows, column bounds retained.
    if (!e.inconsistent_bound_columns.empty()) {
        e.relaxation_status = "not_applicable_column_bounds_inconsistent";
    } else {
        Model r = m;
        r.maximize = false;
        r.obj_const = 0;
        std::fill(r.cost.begin(), r.cost.end(), 0.0);
        r.is_int.assign(ncols, 0);
        std::vector<int> lower_col(nrows, -1), upper_col(nrows, -1);
        std::vector<double> weight(nrows, 1.0);
        for (size_t i = 0; i < nrows; ++i) {
            double mag = 0;
            if (std::isfinite(m.row_lo[i])) mag = std::max(mag, std::abs(m.row_lo[i]));
            if (std::isfinite(m.row_up[i])) mag = std::max(mag, std::abs(m.row_up[i]));
            weight[i] = 1.0 / (1.0 + mag);
            if (std::isfinite(m.row_lo[i])) {  // A x + sl >= lo : the lower side relaxed by sl
                lower_col[i] = int(r.cols.size());
                r.col_names.push_back("RELAX_LO_" + m.row_names[i]);
                r.cols.push_back({{int(i), 1.0}});
                r.cost.push_back(weight[i]);
                r.col_lo.push_back(0); r.col_up.push_back(kInf);
            }
            if (std::isfinite(m.row_up[i])) {  // A x - su <= up : the upper side relaxed by su
                upper_col[i] = int(r.cols.size());
                r.col_names.push_back("RELAX_UP_" + m.row_names[i]);
                r.cols.push_back({{int(i), -1.0}});
                r.cost.push_back(weight[i]);
                r.col_lo.push_back(0); r.col_up.push_back(kInf);
            }
        }
        r.is_int.assign(r.cols.size(), 0);
        ++solves;
        Result rr = solve_lp_gated(r, r.col_lo, r.col_up, nullptr, std::max(left(), 1.0), false, true);
        if (rr.status == Status::Optimal && rr.x.size() == r.cols.size()) {
            e.relaxation_status = "optimal";
            e.relaxation_objective = rr.objective;
            e.relaxation_gap = rr.gap;
            e.relaxation_quality = rr.certificate_quality;
            e.relaxation_x.assign(rr.x.begin(), rr.x.begin() + ncols);
            e.relax_lower.assign(nrows, 0.0);
            e.relax_upper.assign(nrows, 0.0);
            e.relax_weight = weight;
            // Independent recheck on the model with each relaxed side moved by its amount.
            std::vector<long double> act(nrows, 0);
            double worst = 0;
            for (size_t j = 0; j < ncols; ++j) {
                if (std::isfinite(m.col_lo[j])) worst = std::max(worst, (m.col_lo[j] - e.relaxation_x[j]) / (1 + std::abs(m.col_lo[j])));
                if (std::isfinite(m.col_up[j])) worst = std::max(worst, (e.relaxation_x[j] - m.col_up[j]) / (1 + std::abs(m.col_up[j])));
                for (const Entry& en : m.cols[j]) act[en.index] += (long double)en.value * e.relaxation_x[j];
            }
            for (size_t i = 0; i < nrows; ++i) {
                if (lower_col[i] >= 0) e.relax_lower[i] = std::max(0.0, rr.x[lower_col[i]]);
                if (upper_col[i] >= 0) e.relax_upper[i] = std::max(0.0, rr.x[upper_col[i]]);
                double lo = m.row_lo[i] - e.relax_lower[i], up = m.row_up[i] + e.relax_upper[i];
                if (std::isfinite(lo)) worst = std::max(worst, (double)((lo - act[i]) / (1 + std::abs(lo))));
                if (std::isfinite(up)) worst = std::max(worst, (double)((act[i] - up) / (1 + std::abs(up))));
            }
            e.relaxation_violation = worst;
            if (!(worst <= kWitnessTol)) e.relaxation_status = "unverified_point";
        } else {
            e.relaxation_status = std::string("not_solved_") + status_name(rr.status);
        }
    }
    if (!e.unproven_rows.empty() || out_of_time)
        finish("reduced_unproven", "the subsystem is infeasible (verified) but " + std::to_string(e.unproven_rows.size()) +
                                       " row(s) lack a verified feasibility witness for their removal");
    else if (e.rows.empty())
        finish("bounds_only", "the column bounds alone are infeasible; no row is needed");
    else
        finish("irreducible", "row-irreducible under retained column bounds");
    return e;
}
