// Build with src/*.cpp except main.cpp. Deterministic budget tests for explain_infeasibility using the
// simulated-elapsed test hook (no dependence on machine speed).
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>
#include "../src/taral.hpp"

static bool exponent_ok(double s) { int ex = 0; return std::frexp(s, &ex) == 0.5; }  // an exact power of two

static Model conflict() {  // x + y >= 10 and x + y <= 4, 0 <= x,y <= 100, plus an irrelevant row
    Model m;
    m.row_names = {"LOW", "HIGH", "SIDE"};
    m.row_lo = {10, -kInf, -kInf};
    m.row_up = {kInf, 4, 50};
    m.col_names = {"X", "Y"};
    m.cols = {{{0, 1}, {1, 1}, {2, 1}}, {{0, 1}, {1, 1}}};
    m.cost = {1, 1};
    m.col_lo = {0, 0};
    m.col_up = {100, 100};
    m.is_int = {0, 0};
    return m;
}

int main() {
    Model m = conflict();
    InfeasibilityExplanation e = explain_infeasibility(m, 10.0);
    assert(e.status == "irreducible" && e.rows.size() == 2 && e.relaxation_status == "optimal");
    // time is gone once the deletion phase ends: the elastic solve must not start, status stays honest
    e = explain_infeasibility(m, 10.0, 1000.0);
    assert(e.relaxation_status == "not_run_time_limit");
    assert(e.relaxation_x.empty());
    assert(e.relaxation_scale == 0 && e.relaxation_quality.empty());  // the scaled LP was never attempted: nothing is claimed about it
    // the elastic solve runs on costs scaled by an exact power of two in [2^20, 2^30]; the reported objective is in original units
    InfeasibilityExplanation f = explain_infeasibility(m, 10.0);
    assert(f.relaxation_scale >= 1048576.0 && f.relaxation_scale <= 1073741824.0 && exponent_ok(f.relaxation_scale));
    assert(std::abs(f.relaxation_objective - 6.0 / 11.0) < 1e-9);  // original units: 6 * min(1/11, 1/5)
    assert(!f.relaxation_quality.empty());
    long with_elastic = explain_infeasibility(m, 10.0).lp_solves;
    assert(e.lp_solves == with_elastic - 1);
    assert(e.status == "irreducible");  // the row proof itself is complete; only the relaxation report is missing
    // a finite row bound near DBL_MAX gives a denormal weight, so 1/wmin overflows to inf: the exponent must clamp to 30 (no int cast of inf)
    Model big = conflict();
    big.row_up[2] = 1.7e308;
    InfeasibilityExplanation g = explain_infeasibility(big, 10.0);
    assert(g.status == "irreducible" && g.relaxation_scale == 1073741824.0);
    // budget already exhausted at the start: no proof, no solve claims
    e = explain_infeasibility(m, 0.0);
    assert(e.status == "no_verified_proof");
    // relaxation_quality_label on hand-built structs (no solver involved)
    InfeasibilityExplanation h;
    h.relaxation_status = "not_solved_numerical_failure"; h.relaxation_scale = 1048576.0;
    assert(std::string(relaxation_quality_label(h)) == "not_available" && h.relaxation_scale == 1048576.0);
    h = InfeasibilityExplanation(); h.relaxation_status = "not_run_time_limit"; h.relaxation_scale = 0;
    assert(std::string(relaxation_quality_label(h)) == "not_run");
    h = InfeasibilityExplanation(); h.relaxation_status = "unverified_point"; h.relaxation_quality = "pass"; h.relaxation_scale = 1048576.0;
    assert(std::string(relaxation_quality_label(h)) == "pass");
    h = InfeasibilityExplanation(); h.relaxation_status = "not_applicable_column_bounds_inconsistent";
    assert(std::string(relaxation_quality_label(h)) == "not_run");
    std::printf("api tests ok\n");
    return 0;
}
