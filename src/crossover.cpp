// Basis from an approximate primal point (the PDHG-to-simplex handoff). Variables are the structurals
// then the row logicals (activity of the row), as in the simplex warm-start vector. The m variables
// farthest from their bounds (relative to 1 + |value|) become basic, provided they are farther than
// `tol` (`interior` reports how many are); every other variable is nonbasic at its nearest finite bound (free ones at zero). When fewer
// than m variables qualify the remaining slots go to the next-farthest ones, and the simplex's own
// singular-basis repair handles whatever that leaves rank deficient. With row multipliers y the leftover slots
// go to the variables with the smallest reduced cost |c_j - a_j'y| (basic-like), not to arbitrary ones; the sign
// of y is taken from the interior variables, whose reduced cost should vanish.
#include <algorithm>
#include <cmath>
#include <numeric>

#include "taral.hpp"

std::vector<char> basis_from_point(const Model& md, const std::vector<double>& x, const std::vector<double>& y,
                                double tol, int* interior) {
    const int n = int(md.col_names.size()), m = int(md.row_names.size()), N = n + m;
    std::vector<double> v(N, 0.0), lo(N), up(N), dist(N);
    for (int j = 0; j < n; ++j) v[j] = x[j], lo[j] = md.col_lo[j], up[j] = md.col_up[j];
    for (int j = 0; j < n; ++j)
        for (const Entry& e : md.cols[j]) v[n + e.index] += e.value * x[j];
    for (int i = 0; i < m; ++i) lo[n + i] = md.row_lo[i], up[n + i] = md.row_up[i];
    for (int j = 0; j < N; ++j) {
        double d = std::min(v[j] - lo[j], up[j] - v[j]);  // +inf for a free variable
        dist[j] = d / (1 + std::abs(v[j]));
    }
    std::vector<double> score = dist;  // larger = more basic-like
    if (!y.empty()) {
        std::vector<double> dp(N, 0.0), dm(N, 0.0);  // reduced costs for +y and -y
        double sp = 0, sm = 0;
        for (int j = 0; j < n; ++j) {
            double a = 0;
            for (const Entry& e : md.cols[j]) a += e.value * y[e.index];
            dp[j] = md.cost[j] - a, dm[j] = md.cost[j] + a;
        }
        for (int i = 0; i < m; ++i) dp[n + i] = -y[i], dm[n + i] = y[i];  // logical is -e_i with zero cost
        for (int j = 0; j < N; ++j)
            if (dist[j] > tol) sp += std::abs(dp[j]), sm += std::abs(dm[j]);
        const std::vector<double>& d = sp <= sm ? dp : dm;
        for (int j = 0; j < N; ++j)
            if (dist[j] <= tol) score[j] = -std::abs(d[j]) / (1 + (j < n ? std::abs(md.cost[j]) : 0.0));
    }
    std::vector<int> order(N);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return score[a] > score[b]; });
    std::vector<char> basis(N, 1);
    for (int j = 0; j < N; ++j) {
        bool lf = std::isfinite(lo[j]), uf = std::isfinite(up[j]);
        if (!lf && !uf) basis[j] = 3;
        else if (!lf) basis[j] = 2;
        else if (!uf) basis[j] = 1;
        else basis[j] = (v[j] - lo[j] <= up[j] - v[j]) ? 1 : 2;
    }
    for (int k = 0; k < m; ++k) basis[order[k]] = 0;
    if (interior) *interior = int(std::count_if(dist.begin(), dist.end(), [&](double d) { return d > tol; }));
    return basis;
}
