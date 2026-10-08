#include <cassert>
#include <cmath>
#include <iostream>
#include "miqp.hpp"
#include "ipm.hpp"
static int calls;
static int mode;
IpmResult ipm_solve(const Model&, const IpmOptions&) {
    IpmResult r;
    ++calls;
    r.status = IpmStatus::Optimal;
    r.x = {0.5}; r.objective = -1;
    if (calls == 1 && mode == 4) { r.x.clear(); return r; }
    if (calls == 1 && mode == 3) { r.status = IpmStatus::NumericalFailure; return r; }
    if (calls == 2 || calls == 3) { r.x = {0}; r.objective = -1; }
    if (calls >= 4) {
        if (mode == 0 || mode == 5) { r.status = IpmStatus::Infeasible; }
        if (mode == 1) { r.status = IpmStatus::NumericalFailure; }
        if (mode == 2) { r.status = IpmStatus::TimeLimit; }
    }
    if (mode == 5 && calls == 3) r.status = IpmStatus::NumericalFailure;
    return r;
}
int main() {
    Model md; md.cols.resize(1); md.cost={0}; md.col_lo={0}; md.col_up={1}; md.is_int={1};
    for (mode=0; mode<=5; ++mode) {
        calls=0;
        auto r=miqp_proto::solve_miqp(md, 10, 100);
        std::cerr << mode << " " << r.status << " " << r.unresolved << " " << r.has_solution << " " << r.nodes << "\n";
        if (mode == 0) { assert(r.status=="optimal" && r.has_solution && !r.unresolved); }
        if (mode == 1 || mode == 2 || mode == 5) {
            assert(r.has_solution && r.unresolved>0);
            assert(r.status == (mode==2 ? "time_limit" : "numerical_failure"));
            assert(r.best_bound <= -1);
        }
        if (mode == 3 || mode == 4) {
            assert(r.status=="numerical_failure" && !r.has_solution && r.unresolved==1);
            assert(std::isinf(r.best_bound) && r.best_bound < 0);
        }
    }
    calls=0; auto cap=miqp_proto::solve_miqp(md,10,0); assert(cap.status=="time_limit");
    md.maximize=true; assert(miqp_proto::solve_miqp(md,10,100).status=="unsupported");
    std::cout << "8 unresolved-node/status regressions passed\n";
}
