// Opt-in supernodal multifrontal LDL^T (TARAL_QP_SUPER=1) for the interior-point KKT system. Same pivot rule as the
// scalar up-looking Ldl in ipm.cpp (sign test against the largest term that formed the pivot, same replacement value),
// blocked dense kernels. Floating-point contraction is off here as in ipm.cpp, so the result does not depend on the
// host's FMA support. Used only when the knob is set and the factor is large enough to pay.
#pragma once
#include <chrono>
#include <vector>
namespace taral_qp_super {
using Clock = std::chrono::steady_clock;
typedef std::vector<double> Vec;
struct SuperImpl;
class SuperLdl {
  public:
    SuperLdl();
    ~SuperLdl();
    // Lower-by-column entry list as in the scalar Ldl (column k holds rows i <= k, diagonal included), already in
    // elimination order, with the elimination tree parent[] and column counts cnt[] (L column counts).
    void analyze(int n, const std::vector<int>& Kp, const std::vector<int>& Ki, const std::vector<signed char>& sign,
                 const std::vector<int>& parent, const std::vector<int>& cnt);
    // V[p] is the value of entry p of Ki. Returns false when the deadline passes.
    bool factor(const Vec& V, Clock::time_point deadline, double piv_tol, double piv_rep, double piv_abs);
    void solve(Vec& x) const; // in the elimination order
    int nreg() const;
    bool all_finite() const; // D and L finite after factor()
  private:
    SuperImpl* impl_;
};
// Post-order of the elimination tree (children before parents, subtrees contiguous). post[new] = old.
std::vector<int> etree_postorder(const std::vector<int>& parent);
}
