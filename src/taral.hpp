// TARAL-LP engine: shared types and module interfaces. C++17 standard library only.
#pragma once
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
    // cols[p] = entries (row, value) of basis column p. Returns false if rank deficient;
    // the unpivoted positions and rows are then left in bad_pos / bad_rows.
    bool factor(int m, const std::vector<std::vector<Entry>>& cols);
    void ftran(std::vector<double>& rhs_rows, std::vector<double>& out_pos) const;  // B x = rhs
    void btran(std::vector<double>& rhs_pos, std::vector<double>& out_rows) const;  // B'y = rhs
    std::vector<int> bad_pos, bad_rows;

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
    double objective = 0;   // includes obj_const
    long iterations = 0;
    std::string message;
};

Result solve(const Model& model, double time_limit_s);
const char* status_name(Status s);
