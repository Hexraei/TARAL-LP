// Build with src/*.cpp except main.cpp. Usage: det_api MODEL.mps
// With a work cap set, an explicit (tiny) wall limit must still stop the solve with time_limit.
#include <cassert>
#include <cstdio>
#include "../src/taral.hpp"

int main(int argc, char** argv) {
    assert(argc > 1);
    Model m = read_mps(argv[1]);
    work_budget().cap = 1000000;
    work_budget().used = 0;
    Result r = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 1e-9, false, true);
    assert(r.status == Status::TimeLimit);
    assert(work_budget().used < 1000);  // stopped by the clock long before the work cap
    work_budget().used = 0;
    Result q = solve_lp_gated(m, m.col_lo, m.col_up, nullptr, 1e9, false, true);  // no wall pressure: work decides
    assert(q.status == Status::Optimal && work_budget().used == q.iterations);
    std::printf("det api tests ok\n");
    return 0;
}
