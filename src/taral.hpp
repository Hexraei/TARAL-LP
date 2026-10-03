// TARAL-LP engine: shared types and module interfaces. C++17 standard library only.
#pragma once
#include <array>
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

// Thrown when the optional deadline passes while reading. Not a ParseError: the free-format to fixed-column
// fallback in read_mps must not retry after a timeout.
struct ParseTimeLimit : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Throws ParseError, or ParseTimeLimit when `deadline` passes while reading (checked every few thousand lines).
Model read_mps(const std::string& path,
               std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max());

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
    long prop_tightened = 0;                   // propagation: integer bounds tightened
    long prop_crossed = 0, prop_crossed_lp_infeasible = 0, prop_pruned = 0;  // nodes whose integer domain came out
                                               // empty; of those the LP confirmed; pruned without the LP
    // --audit-prop: nodes with an emptied domain, re-checked (see Search::audit_pruned in milp.cpp):
    // {nodes, LP infeasible, LP feasible, LP other, integer-empty by plain B&B, integer-feasible (a bug),
    //  undecided, justified by the incumbent cutoff only, reduced-cost fixings checked, fixings found wrong}
    std::array<long, 10> audit{};
    long rc_fixed = 0, rc_skipped = 0;  // reduced-cost fixing: bounds tightened, bases not trusted
    std::string message;
};
// Warm-start basis from an approximate primal point (src/crossover.cpp); *interior = variables farther than tol from a bound.
std::vector<char> basis_from_point(const Model& model, const std::vector<double>& x,
                                const std::vector<double>& row_duals, double tol, int* interior);
struct MilpOptions {
    bool audit_prop = false;  // --audit-prop: re-check every node whose propagated integer domain came out empty
    bool no_prop_prune = false;  // --no-prop-prune: do not prune such a node directly; the LP decides
};
MilpResult solve_milp(const Model& model, double time_limit_s, long node_limit, MilpOptions opt = {});
const char* status_name(Status s);

// Recomputes proof validity from model coefficients, not solver basis/pricing state.
// Invalid/missing/nonfinite proofs become NumericalFailure, never a proved status.
bool verify_nonoptimal_certificate(const Model&, const std::vector<double>&,
                                  const std::vector<double>&, Result&);
const char* status_name(Status s);

// KKT-gated LP solve (src/kkt_gate.cpp): Optimal only when an independent check on the original model agrees,
// with one power-of-two Ruiz-equilibrated retry before numerical_failure. use_dual selects solve_lp_dual.
struct KktReport {
    bool ok = false;
    double primal = 0, dual = 0;  // worst primal violation, worst dual violation relative to its column cost scale
    std::string msg;
};
KktReport kkt_check(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                    const Result& r);
Result solve_lp_gated(const Model& model, const std::vector<double>& col_lo, const std::vector<double>& col_up,
                      const std::vector<char>* warm_basis, double time_limit_s, bool use_dual, bool primal_fallback = true);

// Verified LP infeasibility explanation (src/explain.cpp, docs/infeasibility_explanation.md). Integrality is
// ignored: the explanation is for the LP relaxation. "Irreducible" means row-irreducible with ALL column
// bounds retained and tolerance-feasible witnesses; it is not a minimum-cardinality set.
struct InfeasibilityExplanation {
    // relaxation_feasible | irreducible | reduced_unproven | bounds_only | no_verified_proof
    std::string status, message;
    std::vector<int> rows;                       // original row indices of the explained subsystem
    std::vector<std::vector<double>> witness;    // per explained row: point feasible for the others (empty if unproven)
    std::vector<double> witness_violation;       // worst relative violation of that witness
    std::vector<int> unproven_rows;
    std::vector<int> bound_columns;              // columns whose bounds carry a nonzero multiplier in the proof
    std::vector<int> inconsistent_bound_columns; // col_lo > col_up
    std::vector<double> farkas_row_lower, farkas_row_upper, farkas_col_lower, farkas_col_upper;  // full model
    bool certificate_verified = false;
    double certificate_margin = 0, certificate_residual = 0;
    // Minimum weighted L1 side relaxation of ALL rows, column bounds retained.
    std::string relaxation_status = "not_run", relaxation_quality;
    double relaxation_objective = 0, relaxation_gap = 0, relaxation_violation = 0;
    std::vector<double> relaxation_x, relax_lower, relax_upper, relax_weight;
    long lp_solves = 0;
    double wall_s = 0;
};
// The budget covers every phase (root test, deletion filter, elastic solve). The last argument is a test hook
// that adds simulated elapsed time after the deletion phase; production callers leave it at 0.
InfeasibilityExplanation explain_infeasibility(const Model& model, double time_limit_s, double test_elapsed_after_deletion_s = 0);

// ---- Verified presolve (src/presolve.cpp) -------------------------------------------------------------
// Opt-in (--presolve). Linear models only (no quadratic objective). Every reduction is recorded in a log
// that an independent replayer (benchmarks/presolve_replay_check.py) re-derives from the ORIGINAL model, and
// every returned point is audited in original space (rows, bounds, integrality, objective) before it is
// reported. Reductions: integer bound rounding, fixed-column substitution, empty/singleton rows,
// activity-redundant rows, empty-column fixing. Floating point throughout, tolerances stated in
// docs/verified_presolve.md.
struct PresolveOp {
    std::string type;       // round_int_bounds, fix_col, empty_row, singleton_row, redundant_row, fix_empty_col
    int row = -1, col = -1;
    double a = 0, v1 = 0, v2 = 0, v3 = 0, v4 = 0;  // type-specific data, see presolve.cpp
};
struct PresolveResult {
    bool infeasible = false;
    std::string infeasible_reason;
    Model reduced;
    std::vector<int> kept_rows, kept_cols;   // original indices of the reduced model's rows/columns
    std::vector<double> fixed_value;         // per original column; valid where the column was removed
    std::vector<char> removed_col;
    std::vector<PresolveOp> log;
    int passes = 0;
    bool timed_out = false;  // deadline hit: the partial reduction is NOT used; the log is kept for inspection
};
// The deadline is checked at the start of every pass and every 256 rows/columns. test_timeout_after_ops >= 0 is a
// test hook that behaves as if the deadline passed once that many reductions are logged.
PresolveResult presolve_model(const Model& orig,
                              std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max(),
                              long test_timeout_after_ops = -1);
std::vector<double> presolve_expand(const Model& orig, const PresolveResult& p, const std::vector<double>& x_reduced);
struct PresolveAudit {
    bool ok = false;
    double max_row_violation = 0, max_bound_violation = 0, max_int_violation = 0;  // relative: viol/(1+|bound|)
    double objective = 0, reported_objective = 0, objective_diff = 0;
};
PresolveAudit presolve_audit(const Model& orig, const std::vector<double>& x, double reported_objective);
bool write_presolve_log(const char* path, const Model& orig, const PresolveResult& p, const PresolveAudit* audit);
