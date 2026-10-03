// Build with src/*.cpp except main.cpp. Deadline discipline of presolve_model (deterministic test hook).
#include <cassert>
#include <chrono>
#include <cstdio>
#include "../src/taral.hpp"

int main() {
    Model m;  // x0 fixed at 1, x1 free in [0,4] with a singleton row, so there are several reductions
    m.row_names = {"R0", "R1"};
    m.row_lo = {-kInf, 1};
    m.row_up = {10, kInf};
    m.col_names = {"X0", "X1"};
    m.cols = {{{0, 1}}, {{0, 1}, {1, 1}}};
    m.cost = {1, 1};
    m.col_lo = {1, 0};
    m.col_up = {1, 4};
    m.is_int = {0, 0};
    PresolveResult p = presolve_model(m);
    assert(!p.timed_out && !p.log.empty());
    // hook: pretend the deadline passes after the first logged reduction
    PresolveResult q = presolve_model(m, std::chrono::steady_clock::time_point::max(), 1);
    assert(q.timed_out && q.log.size() == 1);
    // already expired deadline: no reduction is logged, state is reported as timed out
    PresolveResult z = presolve_model(m, std::chrono::steady_clock::now() - std::chrono::seconds(1));
    assert(z.timed_out && z.log.empty());
    std::printf("presolve api tests ok\n");
    return 0;
}
