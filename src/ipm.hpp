// Primal-dual interior-point method for LP and convex QP with linear constraints.
// C++17 standard library only. See ipm.cpp for the algorithm and its limits.
#pragma once
#include <string>
#include <vector>

#include "taral.hpp"

enum class IpmStatus {
    Optimal,        // KKT measures below tol on the original model (re-checked from scratch)
    Infeasible,     // approximate Farkas certificate verified, or a trivially infeasible row/bound
    Unbounded,      // recession direction verified and primal residual already small
    DualInfeasible, // recession direction verified, primal feasibility not established
    TimeLimit,
    IterationLimit,
    NumericalFailure, // stalled, or converged internally but failed the original-model recheck
    Nonconvex,        // Q (in minimisation sense) has a negative pivot: not positive semidefinite
    Unsupported,      // integer columns
};

struct IpmOptions {
    double tol = 1e-8;      // relative primal residual, dual residual and gap on the original model
    double time_limit = 60; // seconds, wall clock
    int max_iter = 200;
    bool verbose = false; // per-iteration log on stderr
};

// Sign convention, in the model's own sense (min or max):  reduced_cost = c + Q x - A' row_dual.
// Measures (all computed on the original model from x and row_dual alone):
//   primal_res = max over rows/columns of violation / (1 + |violated bound|)
//   dual_res   = max sign violation of row_dual / reduced_cost w.r.t. which bounds are finite, / (1 + ||c||inf)
//   gap        = |primal obj - dual obj| / (1 + |primal obj|), dual obj = bound terms - 0.5 x'Qx + obj_const
struct IpmResult {
    IpmStatus status = IpmStatus::NumericalFailure;
    std::vector<double> x, row_activity, row_dual, reduced_cost;
    double objective = 0, dual_objective = 0; // include obj_const, model's sense
    long iterations = 0;
    double primal_res = 0, dual_res = 0, gap = 0;
    double max_row_viol = 0, max_bound_viol = 0; // absolute
    long long factor_nnz = 0;                    // nonzeros in L of the augmented system
    int kkt_dim = 0;
    std::string message;
};

IpmResult ipm_solve(const Model& model, const IpmOptions& options = IpmOptions());
const char* ipm_status_name(IpmStatus s);
