// Unit tests for src/kkt_gate.cpp. Build: g++ -std=c++17 -O2 -I../src kkt_gate_tests.cpp ../src/{mps,simplex,dual,lu,kkt_gate}.cpp -o kkt_gate_tests
#include <cassert>
#include <cmath>
#include <cstdio>
#include "taral.hpp"

static Model two_var() {  // min x + 2y  s.t. x + y >= 2, x,y >= 0.  optimum (2,0), objective 2
    Model m; m.col_names = {"x", "y"}; m.row_names = {"r"};
    m.cols = {{{0, 1}}, {{0, 1}}}; m.cost = {1, 2}; m.col_lo = {0, 0}; m.col_up = {kInf, kInf};
    m.row_lo = {2}; m.row_up = {kInf};
    return m;
}
int main() {
    int passed = 0;
    { // 1. a true optimum passes the gate
        Model m = two_var(); Result r = solve_lp(m, m.col_lo, m.col_up, nullptr, 10);
        assert(r.status == Status::Optimal); KktReport k = kkt_check(m, m.col_lo, m.col_up, r);
        assert(k.ok); Result g = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 10, false);
        assert(g.status == Status::Optimal && std::abs(g.objective - 2) < 1e-9); ++passed; }
    { // 2. a feasible but non-optimal point carrying the optimal basis flags is rejected (dual side, not primal)
        Model m = two_var(); Result r = solve_lp(m, m.col_lo, m.col_up, nullptr, 10);
        Result bad = r; bad.x = {0, 2}; bad.objective = 4;  // feasible: 0+2 >= 2
        KktReport k = kkt_check(m, m.col_lo, m.col_up, bad); assert(!k.ok); ++passed; }
    { // 3. the "optimal at 0" failure mode: an infeasible point claimed optimal is rejected
        Model m = two_var(); Result r = solve_lp(m, m.col_lo, m.col_up, nullptr, 10);
        Result bad = r; bad.x = {0, 0}; bad.objective = 0;
        KktReport k = kkt_check(m, m.col_lo, m.col_up, bad); assert(!k.ok); ++passed; }
    { // 4. a primal-feasible point whose basis gives wrong-sign reduced costs is rejected (cost scale judgement)
        Model m = two_var(); Result r = solve_lp(m, m.col_lo, m.col_up, nullptr, 10);
        Result bad = r; bad.x = {0, 2}; bad.basis = r.basis; std::swap(bad.basis[0], bad.basis[1]);
        KktReport k = kkt_check(m, m.col_lo, m.col_up, bad); assert(!k.ok); ++passed; }
    { // 5. missing basis cannot be certified
        Model m = two_var(); Result r = solve_lp(m, m.col_lo, m.col_up, nullptr, 10); r.basis.clear();
        assert(!kkt_check(m, m.col_lo, m.col_up, r).ok); ++passed; }
    { // 6. badly scaled copy of the same LP (columns x1e7 and x1e-7, rows x1e6): gated result matches the known optimum
        Model m = two_var(); double sx = 1e7, sy = 1e-7, rr = 1e6;
        for (Entry& e : m.cols[0]) e.value *= rr * sx;
        for (Entry& e : m.cols[1]) e.value *= rr * sy;
        m.cost[0] *= sx; m.cost[1] *= sy; m.row_lo[0] *= rr;  // x' = x/sx, y' = y/sy
        Result g = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 10, false);
        assert(g.status == Status::Optimal && std::abs(g.objective - 2) < 1e-6 * 2);
        Result d = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 10, true);
        assert(d.status == Status::Optimal && std::abs(d.objective - 2) < 1e-6 * 2); ++passed; }
    { // 7. maximize sense
        Model m = two_var(); m.maximize = true; m.cost = {-1, -2};  // max -x-2y == min x+2y, objective -2
        Result g = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 10, false);
        assert(g.status == Status::Optimal && std::abs(g.objective + 2) < 1e-9);
        assert(kkt_check(m, m.col_lo, m.col_up, g).ok); ++passed; }
    { // 8. infeasible and unbounded statuses pass through unchanged
        Model m = two_var(); m.row_up = {1}; m.row_lo = {2};
        assert(solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 10, false).status == Status::Infeasible);
        Model u = two_var(); u.cost = {-1, -2}; assert(solve_lp_gated(u, u.col_lo, u.col_up, nullptr, 10, false).status == Status::Unbounded); ++passed; }
    std::printf("kkt_gate_tests: %d/8 passed\n", passed);
    return passed == 8 ? 0 : 1;
}
