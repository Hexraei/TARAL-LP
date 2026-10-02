// Build separately with src/{mps,simplex,dual,lu,nonoptimal_certificate}.cpp.
#include <cassert>
#include <cmath>
#include "../src/taral.hpp"
int main() {
    Model m; m.col_names={"x"}; m.cols={{}}; m.cost={-1}; m.col_lo={0}; m.col_up={kInf};
    Result r=solve(m,10); assert(r.status==Status::Unbounded && r.certificate_verified);
    r.ray[0]=-1; assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    assert(r.status==Status::NumericalFailure && !r.certificate_verified);
    r=solve(m,10);r.x[0]=NAN; assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    r=solve(m,10);r.ray.clear(); assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    m.col_lo={2};m.col_up={1};r=solve(m,10);assert(r.status==Status::Infeasible);
    r.farkas_col_lower[0]=0;assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    r=solve(m,10);r.farkas_col_upper[0]=-1;assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    r=solve(m,10);m.cost[0]=NAN;assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    // Exported verification must work against the ORIGINAL maximize model,
    // not just the internally negated minimization problem.
    m.cost={1};m.col_lo={0};m.col_up={kInf};m.maximize=true;
    r=solve(m,10);assert(r.status==Status::Unbounded && r.certificate_verified);
    assert(verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r) && r.certificate_margin==0.5);
    r=solve_lp_dual(m,m.col_lo,m.col_up,nullptr,10);
    assert(r.status==Status::Unbounded && verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    r.ray[0]=-1;assert(!verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    m.cost={-1};m.col_lo={-kInf};m.col_up={0};r=solve(m,10);
    assert(r.status==Status::Unbounded && verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    m.col_lo={2};m.col_up={1};r=solve(m,10);
    assert(r.status==Status::Infeasible && verify_nonoptimal_certificate(m,m.col_lo,m.col_up,r));
    m.maximize=false;
    // Near-cancelled Farkas multipliers on a directly free column: original
    // singleton rows imply finite endpoints that safely bound roundoff.
    Model box;box.col_names={"free"};box.cols={{{0,1},{1,1}}};box.cost={0};
    box.col_lo={-kInf};box.col_up={kInf};box.row_names={"lower","upper"};
    box.row_lo={2,-kInf};box.row_up={kInf,1};
    Result proof;proof.status=Status::Infeasible;
    proof.farkas_row_lower={0.5+1e-14,0};proof.farkas_row_upper={0,0.5-1e-14};
    proof.farkas_col_lower={0};proof.farkas_col_upper={0};
    assert(verify_nonoptimal_certificate(box,box.col_lo,box.col_up,proof));
    // If singleton evidence is removed, the residual is not bounded: reject.
    box.col_names.push_back("other");box.cols.push_back({{0,1},{1,1}});box.cost.push_back(0);
    box.col_lo.push_back(-kInf);box.col_up.push_back(kInf);
    proof.status=Status::Infeasible;proof.farkas_col_lower.push_back(0);proof.farkas_col_upper.push_back(0);
    assert(!verify_nonoptimal_certificate(box,box.col_lo,box.col_up,proof));
    // Duplicate (row,column) entries must be aggregated before bound inference:
    // this FEASIBLE model (x=1e8) must not accept a fake Farkas proof.
    Model dup;dup.col_names={"x"};dup.cols={{{0,1},{1,2},{1,-1}}};dup.cost={0};
    dup.col_lo={-kInf};dup.col_up={kInf};dup.row_names={"lower","upper","zero"};
    dup.row_lo={1e8,-kInf,0};dup.row_up={kInf,1e8,kInf};
    Result fake;fake.status=Status::Infeasible;fake.farkas_row_lower={1e-14,0,1-1e-14};
    fake.farkas_row_upper={0,0,0};fake.farkas_col_lower={0};fake.farkas_col_upper={0};
    assert(!verify_nonoptimal_certificate(dup,dup.col_lo,dup.col_up,fake));
    assert(!fake.certificate_verified && fake.certificate_margin==0);
    // A row with two nonzero aggregated columns is not a singleton even if one entry cancels.
    Model cancel;cancel.col_names={"x","y"};cancel.cols={{{0,1},{1,1},{1,-1}},{{0,1}}};cancel.cost={0,0};
    cancel.col_lo={-kInf,-kInf};cancel.col_up={kInf,kInf};cancel.row_names={"r0","r1"};
    cancel.row_lo={1e8,-kInf};cancel.row_up={kInf,1e8};
    Result fake2;fake2.status=Status::Infeasible;fake2.farkas_row_lower={1e-14,0};
    fake2.farkas_row_upper={0,1-1e-14};fake2.farkas_col_lower={0,0};fake2.farkas_col_upper={0,0};
    assert(!verify_nonoptimal_certificate(cancel,cancel.col_lo,cancel.col_up,fake2));
    // Roundoff residual on a column with no finite bound on the needed side. (a) The
    // column is free, but an equality row ties it to a boxed column, so its feasible
    // values are bounded: the implied bound compensates the residual and the proof is valid.
    Model tie;tie.col_names={"x","y"};tie.cols={{{1,1}},{{0,1},{1,-1}}};tie.cost={0,0};
    tie.col_lo={-kInf,0};tie.col_up={kInf,1};tie.row_names={"r","def"};
    tie.row_lo={2,0};tie.row_up={kInf,0};
    Result tp;tp.status=Status::Infeasible;tp.farkas_row_lower={0.5,1e-14};tp.farkas_row_upper={0,0};
    tp.farkas_col_lower={0,0};tp.farkas_col_upper={0,0.5};
    assert(verify_nonoptimal_certificate(tie,tie.col_lo,tie.col_up,tp));
    // Same residual with the tie row one-sided (x >= y only), so x has no implied upper
    // bound and the residual is not bounded: must be rejected.
    tie.row_up[1]=kInf;tp.status=Status::Infeasible;
    assert(!verify_nonoptimal_certificate(tie,tie.col_lo,tie.col_up,tp));
    // (b) A lower-unbounded column with a covering upper multiplier: the residual is
    // removed exactly by lowering that multiplier, so the proof stays valid.
    Model up1;up1.col_names={"x","w"};up1.cols={{{0,1}},{{0,1}}};up1.cost={0,0};
    up1.col_lo={-kInf,0};up1.col_up={0,0};up1.row_names={"r"};up1.row_lo={1};up1.row_up={kInf};
    Result ur;ur.status=Status::Infeasible;ur.farkas_row_lower={1.0/3-1e-14};ur.farkas_row_upper={0};
    ur.farkas_col_lower={0,0};ur.farkas_col_upper={1.0/3,1.0/3-1e-14};
    assert(verify_nonoptimal_certificate(up1,up1.col_lo,up1.col_up,ur));
    // The same residual WITHOUT a covering multiplier is not repairable.
    ur.status=Status::Infeasible;ur.farkas_col_upper[0]=0;
    assert(!verify_nonoptimal_certificate(up1,up1.col_lo,up1.col_up,ur));
    // Budget exhaustion is not a proof.
    m.cost={-1};m.col_lo={0};m.col_up={kInf};assert(solve(m,-1).status==Status::TimeLimit);
    // Too-small strict margin is explicitly unknown/numerical failure.
    m.cost={-1e-9};r=solve(m,10);assert(r.status!=Status::Unbounded);
}
