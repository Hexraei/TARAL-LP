#include <cassert>
#include <cmath>
#include "../src/taral.hpp"
int main() {
    Model m;m.row_names={"zero"};m.row_lo={0};m.row_up={0};m.col_names={"x"};
    m.cols={{{0,1}}};m.cost={0};m.col_lo={0};m.col_up={0};
    Result r=solve(m,10);
    assert(r.status==Status::Optimal && r.certificate_quality=="pass");
    assert(r.row_term_magnitude[0]==0 && r.row_violation_magnitude_scaled[0]==0);
    assert(solve(m,-1).certificate_quality=="unknown");
    m.row_lo={1};m.row_up={1};m.col_lo={1};m.col_up={1};r=solve(m,10);
    assert(r.status==Status::Optimal && r.certificate_quality=="pass");
    assert(r.row_term_magnitude[0]==2 && r.row_violation_abs[0]==0);
    // Violation smaller than phase-1 tolerance is not an infeasibility proof;
    // check the magnitude ratio directly from original-row raw activity.
    m.row_lo={0};m.row_up={0};m.col_lo={5e-10};m.col_up=m.col_lo;r=solve(m,10);
    assert(r.status==Status::Optimal && r.certificate_quality=="pass");
    assert(r.row_violation_abs[0]==5e-10 && r.row_term_magnitude[0]==5e-10);
    assert(r.row_violation_magnitude_scaled[0]==1);
    // Endpoint magnitude is continuous at roundoff: an activity a few ulps below or above a
    // range endpoint reports the same row_term_magnitude (the endpoint counts in both).
    double below=0,above=0;
    for (double x : {1-1e-15, 1+1e-15}) {
        Model t;t.row_names={"r"};t.row_lo={1};t.row_up={9};t.col_names={"x"};
        t.cols={{{0,1}}};t.cost={0};t.col_lo={x};t.col_up={x};
        Result q=solve(t,10);assert(q.status==Status::Optimal);
        (x<1?below:above)=q.row_term_magnitude[0];
        assert(q.row_violation_abs[0]==(x<1?1-x:0));
    }
    assert(std::abs(below-above)<1e-12 && std::abs(below-2)<1e-12);
}
