#pragma once
#include <cmath>
// R12: round a node bound b UP to the objective lattice {k + off : k integer}, weakened for safety. Valid only when cost'x is
// integer-valued for every integer-feasible x (checked from the model by the caller, which also requires |off| <= 1e9).
// Result r <= (exact) ceil(b - off - eps) + off for eps = 1e-6(1+|b|), also under floating-point cancellation: the slack added to
// the subtraction covers its rounding error (<= 1 ulp of max(|b|,|off|)), and r is stepped down by ulps of its own operands.
// This is a numerical-experiment guard, not a proof about the LP bound b itself (see README-R12.md limitations).
inline double objint_round_up(double b, double off) {
    if (!std::isfinite(b)) return b;
    const double m = std::max(std::abs(b), std::abs(off));
    const double ulp_m = std::nextafter(m, INFINITY) - m;
    const double k = std::ceil((b - off) - (1e-6 * (1.0 + std::abs(b)) + 8.0 * ulp_m));
    double r = k + off;
    if (off != 0.0) {
        const double m2 = std::max(std::abs(k), std::abs(off));
        r -= 8.0 * (std::nextafter(m2, INFINITY) - m2);
    }
    return r;
}
