// Build with src/*.cpp except main.cpp. Deterministic budget tests for explain_infeasibility using the
// simulated-elapsed test hook (no dependence on machine speed).
#include <cassert>
#include <cmath>
#include <cstdio>
#include "../src/taral.hpp"

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
    long with_elastic = explain_infeasibility(m, 10.0).lp_solves;
    assert(e.lp_solves == with_elastic - 1);
    assert(e.status == "irreducible");  // the row proof itself is complete; only the relaxation report is missing
    // budget already exhausted at the start: no proof, no solve claims
    e = explain_infeasibility(m, 0.0);
    assert(e.status == "no_verified_proof");
    std::printf("api tests ok\n");
    return 0;
}
