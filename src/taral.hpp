// TARAL-LP engine: shared types and module interfaces. C++17 standard library only.
#pragma once
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

constexpr double kInf = std::numeric_limits<double>::infinity();

struct Entry {
    int index;
    double value;
};

struct QEntry {
    int row, col;  // row >= col; Q is symmetric, so each off-diagonal pair is stored once
    double value;
};

// min (or max) cost'x + 0.5 x'Qx + obj_const  subject to  row_lo <= A x <= row_up,  col_lo <= x <= col_up,
// x_j integer where is_int[j].
struct Model {
    std::string name;
    std::vector<std::string> row_names, col_names;
    std::vector<std::vector<Entry>> cols;  // column-wise A: (row, value), rows ascending
    std::vector<double> cost, col_lo, col_up, row_lo, row_up;
    double obj_const = 0;   // minus the RHS given on the objective row
    bool maximize = false;  // OBJSENSE MAX
    std::vector<char> is_int;   // MARKER INTORG..INTEND columns and BV/LI/UI bounds
    std::vector<QEntry> qobj;   // QUADOBJ / QMATRIX / QSECTION on the objective row
    bool has_integers() const;
};

struct ParseError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

Model read_mps(const std::string& path);  // throws ParseError

// Sparse LU of a square basis matrix with Markowitz pivoting and a column threshold.
class SparseLU {
public:
    using Clock = std::chrono::steady_clock;
    // cols[p] = entries (row, value) of basis column p. Returns false if rank deficient;
    // the unpivoted positions and rows are then left in bad_pos / bad_rows. Returns false with
    // timed_out set (bad_pos empty) when `deadline` passes first; the default never expires.
    bool factor(int m, const std::vector<std::vector<Entry>>& cols, Clock::time_point deadline = Clock::time_point::max());
    void ftran(std::vector<double>& rhs_rows, std::vector<double>& out_pos) const;  // B x = rhs
    void btran(std::vector<double>& rhs_pos, std::vector<double>& out_rows) const;  // B'y = rhs
    std::vector<int> bad_pos, bad_rows;
    bool timed_out = false;
    int pivots_done = 0;       // set on timeout: pivots completed, and entries left in the active submatrix
    long active_nnz = 0;

private:
    int m_ = 0;
    std::vector<int> prow_, pcol_;
    std::vector<double> piv_;
    std::vector<std::vector<Entry>> L_, U_;  // L_[k]: (row, multiplier); U_[k]: (position, value)
};

enum class Status { Optimal, Infeasible, Unbounded, TimeLimit, IterationLimit, NumericalFailure };

struct Result {
    Status status = Status::NumericalFailure;
    std::vector<double> x;  // structural values, valid when Optimal
    std::string certificate_quality = "unknown";  // reporting only, strict 1e-8 KKT
    std::vector<double> row_violation_abs, row_violation_magnitude_scaled, row_term_magnitude;
    double max_row_violation_magnitude_scaled = 0;
    double objective = 0;   // includes obj_const
    // Nonoptimal LP proofs, checked against the input rows and effective column bounds.
    std::vector<double> farkas_row_lower, farkas_row_upper, farkas_col_lower, farkas_col_upper;
    std::vector<double> ray;  // recession direction; x is a verified feasible anchor
    bool certificate_verified = false;
    double certificate_residual = 0, certificate_margin = 0;
    long iterations = 0;
    std::string message;
    // LP optimality certificate in model sense: reduced_cost = cost - A' row_dual.
    // Multipliers in minimisation sense are positive at lower bounds, negative at upper
    // bounds; maximise reverses these signs. solve_lp uses its supplied column bounds.
    std::vector<double> row_activity, row_dual, reduced_cost;
    double dual_objective = 0;  // includes obj_const
    double primal_res = 0, dual_res = 0, gap = 0, complementarity = 0;
    double max_row_viol = 0, max_bound_viol = 0;  // absolute violations
    // Residuals: violation/(1+|violated bound|), sign violation/(1+||cost||inf),
    // gap and max absolute multiplier*slack divided by (1+|objective|).
    std::vector<char> basis;  // when Optimal: per variable (structurals, then row logicals)
                              // 0 basic, 1 at lower, 2 at upper, 3 free at zero
};

Result solve(const Model& model, double time_limit_s);
// solve() with the column bounds replaced by col_lo/col_up and, optionally, a warm start from the
// basis of an earlier solve of the same model (any bounds).
Result solve_lp(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                const std::vector<char>* warm_basis, double time_limit_s);
// Same contract, solved by the bounded dual simplex (src/dual.cpp) followed by a primal clean-up
// from its basis that restores the true costs and certifies the result.
Result solve_lp_dual(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                     const std::vector<char>* warm_basis, double time_limit_s);

// Branch and bound over solve_lp. Objective values are in the model's own sense.
struct MilpResult {
    std::string status;  // optimal infeasible unbounded unbounded_relaxation time_limit node_limit numerical_failure
    bool has_solution = false;
    std::vector<double> x;          // incumbent, checked against the original model
    double objective = 0;           // incumbent objective
    double best_bound = 0;          // proven bound (-inf/+inf when none)
    double gap = kInf;              // |objective - best_bound| / max(1, |objective|)
    long nodes = 0, lp_iterations = 0, unresolved_nodes = 0;
    std::string message;
};
MilpResult solve_milp(const Model& model, double time_limit_s, long node_limit);
// Warm-start basis from an approximate primal point (src/crossover.cpp); *interior = variables farther than tol from a bound.
std::vector<char> basis_from_point(const Model& model, const std::vector<double>& x,
                                const std::vector<double>& row_duals, double tol, int* interior);
const char* status_name(Status s);

// Recomputes proof validity from model coefficients, not solver basis/pricing state.
// Invalid/missing/nonfinite proofs become NumericalFailure, never a proved status.
bool verify_nonoptimal_certificate(const Model&, const std::vector<double>&,
                                  const std::vector<double>&, Result&);
