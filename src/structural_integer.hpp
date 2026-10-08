// Exact, restricted structural integer proofs. No solver library dependencies.
// All arithmetic decisions use exactly represented integral doubles with bounded
// magnitude. Unsupported structures fall through without changing the model.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <vector>
#include "taral.hpp"

namespace structural_integer {
using Bits = std::vector<uint64_t>;
inline bool bit(const Bits& a, int k) { return (a[k / 64] >> (k % 64)) & 1; }
inline void toggle(Bits& a, int k) { a[k / 64] ^= uint64_t(1) << (k % 64); }
inline void add(Bits& a, const Bits& b) {
    for (size_t i = 0; i < a.size(); ++i) a[i] ^= b[i];
}
inline bool exact_int(double x) { return std::isfinite(x) && std::abs(x) <= 1099511627776.0 && x == std::trunc(x); }
struct Proof {
    enum Kind { None, Infeasible, Optimal, Timeout } kind = None;
    std::vector<double> x;
    std::string json;
};
inline std::string rows_json(const Bits& b, const std::vector<int>& eq) {
    std::ostringstream s;
    s << '[';
    bool comma = false;
    for (size_t i = 0; i < eq.size(); ++i)
        if (bit(b, int(i))) {
            if (comma) s << ',';
            s << eq[i];
            comma = true;
        }
    s << ']';
    return s.str();
}
// At most 2048 equalities/variables. A parity contradiction is valid even when
// the continuous model is feasible. A unique binary assignment is NOT by itself
// a feasible point: remaining integer auxiliaries and every original row/bound
// must be checked exactly before declaring optimality.
inline Proof parity(const Model& m, std::chrono::steady_clock::time_point deadline) {
    Proof out;
    const int n = int(m.cols.size()), nr = int(m.row_lo.size());
    if (n == 0 || n > 2048 || nr > 4096 || !m.qobj.empty() || int(m.is_int.size()) != n) return out;
    for (int j = 0; j < n; ++j)
        if (!m.is_int[j]) return out;
    std::vector<std::vector<Entry>> rows(nr);
    for (int j = 0; j < n; ++j)
        for (const auto& e : m.cols[j]) rows[e.index].push_back({j, e.value});
    std::vector<int> eq;
    for (int i = 0; i < nr; ++i)
        if (m.row_lo[i] == m.row_up[i] && exact_int(m.row_lo[i])) {
            bool safe = true;
            for (const auto& e : rows[i]) safe &= exact_int(e.value);
            if (safe) eq.push_back(i);
        }
    if (eq.empty() || eq.size() > 2048) return out;
    const int ne = int(eq.size());
    std::vector<Bits> a(ne, Bits((n + 1 + 63) / 64)), why(ne, Bits((ne + 63) / 64));
    for (int i = 0; i < ne; ++i) {
        for (const auto& e : rows[eq[i]])
            if (int64_t(e.value) % 2) toggle(a[i], e.index);
        if (int64_t(m.row_lo[eq[i]]) % 2) toggle(a[i], n);
        toggle(why[i], i);
    }
    int rank = 0;
    std::vector<int> piv;
    for (int j = 0; j < n; ++j) {
        if (std::chrono::steady_clock::now() >= deadline) {
            out.kind = Proof::Timeout;
            return out;
        }
        int k = rank;
        while (k < ne && !bit(a[k], j)) ++k;
        if (k == ne) continue;
        std::swap(a[k], a[rank]);
        std::swap(why[k], why[rank]);
        for (int i = 0; i < ne; ++i)
            if (i != rank && bit(a[i], j)) {
                add(a[i], a[rank]);
                add(why[i], why[rank]);
            }
        piv.push_back(j);
        ++rank;
    }
    for (int i = rank; i < ne; ++i)
        if (bit(a[i], n)) {
            // Independently recompute XOR from original equations. Never trust only
            // the eliminated matrix when producing an infeasibility status.
            Bits check((n + 1 + 63) / 64);
            for (int k = 0; k < ne; ++k)
                if (bit(why[i], k)) {
                    for (const auto& e : rows[eq[k]])
                        if (int64_t(e.value) % 2) toggle(check, e.index);
                    if (int64_t(m.row_lo[eq[k]]) % 2) toggle(check, n);
                }
            bool ok = bit(check, n);
            for (int j = 0; j < n; ++j) ok &= !bit(check, j);
            if (!ok) return out;
            out.kind = Proof::Infeasible;
            out.json = "{\"type\":\"integer_parity_contradiction\",\"rows\":" + rows_json(why[i], eq) + "}";
            return out;
        }
    // Solve only the fully determined binary + one-row auxiliary pattern.
    // General integers have a residue, not a value, so never fix them from GF(2).
    std::vector<char> binary(n, 0), known(n, 0);
    std::vector<double> x(n, 0);
    for (int j = 0; j < n; ++j) binary[j] = m.col_lo[j] == 0.0 && m.col_up[j] == 1.0;
    std::ostringstream cert;
    cert << "{\"type\":\"integer_parity_unique\",\"fixed_binary\":[";
    bool comma = false;
    for (int i = 0; i < rank; ++i) {
        int j = piv[i];
        if (!binary[j]) continue;
        bool unique = true;
        for (int k = 0; k < n; ++k)
            if (k != j && bit(a[i], k)) unique = false;
        if (!unique) continue;
        x[j] = bit(a[i], n) ? 1.0 : 0.0;
        known[j] = 1;
        if (comma) cert << ',';
        comma = true;
        cert << "{\"col\":" << j << ",\"value\":" << int(x[j]) << ",\"rows\":" << rows_json(why[i], eq) << '}';
    }
    cert << "]}";
    for (int j = 0; j < n; ++j)
        if (binary[j] && !known[j]) return out;
    for (int j = 0; j < n; ++j)
        if (!known[j]) {
            if (m.cols[j].size() != 1) return out;
            auto e = m.cols[j][0];
            int i = e.index;
            if (m.row_lo[i] != m.row_up[i] || !exact_int(m.row_lo[i]) || !exact_int(e.value) || e.value == 0) return out;
            int64_t total = 0;
            for (const auto& t : rows[i])
                if (t.index != j) {
                    if (!known[t.index] || !binary[t.index] || !exact_int(t.value)) return out;
                    // Each value is binary here. 2048 * 2^40 stays safely in int64_t.
                    total += int64_t(t.value) * int64_t(x[t.index]);
                }
            int64_t rhs = int64_t(m.row_lo[i]) - total, c = int64_t(e.value);
            if (rhs % c) return out; // valid necessary congruence but unsupported full solution
            int64_t v = rhs / c;
            if (std::abs(double(v)) > 1099511627776.0) return out;
            x[j] = double(v);
            known[j] = 1;
        }
    // EXACT row/bound validation for this restricted small-integer input.
    for (int j = 0; j < n; ++j) {
        if (!exact_int(x[j]) || x[j] < m.col_lo[j] || x[j] > m.col_up[j]) return out;
        // Restrict auxiliaries further so sum of products fits int64_t.
        if (std::abs(x[j]) > 1048576.0) return out;
    }
    for (int i = 0; i < nr; ++i) {
        int64_t sum = 0;
        for (const auto& e : rows[i]) {
            if (!exact_int(e.value) || std::abs(e.value) > 1048576.0) return out;
            sum += int64_t(e.value) * int64_t(x[e.index]);
        }
        if (std::isfinite(m.row_lo[i]) && (!exact_int(m.row_lo[i]) || double(sum) < m.row_lo[i])) return out;
        if (std::isfinite(m.row_up[i]) && (!exact_int(m.row_up[i]) || double(sum) > m.row_up[i])) return out;
    }
    // All variables are unique: objective coefficients need not be integral.
    out.kind = Proof::Optimal;
    out.x = std::move(x);
    out.json = cert.str();
    return out;
}

inline int64_t gcd_inverse(int64_t a, int64_t mod) {
    int64_t oldr = a, r = mod, oldt = 1, t = 0;
    while (r) {
        int64_t q = oldr / r, nr = oldr - q * r;
        oldr = r;
        r = nr;
        int64_t nt = oldt - q * t;
        oldt = t;
        t = nt;
    }
    if (oldr != 1) return -1;
    oldt %= mod;
    if (oldt < 0) oldt += mod;
    return oldt;
}
// Restricted exact three-integer equality. Enumeration is capped and deadline
// checked; unsupported structures or exhausted enumeration fall back unchanged.
inline Proof three_integer(const Model& m, std::chrono::steady_clock::time_point deadline) {
    Proof out;
    if (m.maximize || m.cols.size() != 3 || m.row_lo.size() != 1 || !m.qobj.empty()) return out;
    if (m.row_lo[0] != 0 || m.row_up[0] != 0 || m.is_int.size() != 3) return out;
    int v = -1;
    for (int j = 0; j < 3; ++j)
        if (m.cost[j] != 0) {
            if (v >= 0 || m.cost[j] != 1) return out;
            v = j;
        }
    if (v < 0 || !exact_int(m.col_lo[v]) || m.col_lo[v] < 1 || m.col_lo[v] > 1000000) return out;
    int vars[2], k = 0;
    int64_t a[3];
    for (int j = 0; j < 3; ++j) {
        if (!m.is_int[j] || m.cols[j].size() != 1 || m.cols[j][0].index != 0) return out;
        double c = m.cols[j][0].value;
        if (!exact_int(c) || std::abs(c) > 1000000 || c == 0) return out;
        a[j] = int64_t(c);
        if (j != v) {
            if (m.col_lo[j] != 0 || m.col_up[j] != kInf) return out;
            vars[k++] = j;
        }
    }
    if (a[v] < 0 || a[vars[0]] > 0 || a[vars[1]] > 0) return out;
    int64_t A = a[v], B = -a[vars[0]], C = -a[vars[1]], inv = gcd_inverse(B, C);
    if (inv < 0 || C <= 1) return out;
    int64_t start = int64_t(m.col_lo[v]), cap = 1000000;
    if (std::isfinite(m.col_up[v])) {
        if (!exact_int(m.col_up[v])) return out;
        cap = std::min(cap, int64_t(m.col_up[v]));
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        out.kind = Proof::Timeout;
        return out;
    }
    for (int64_t x = start; x <= cap; ++x) {
        if ((x & 255) == 0 && std::chrono::steady_clock::now() >= deadline) {
            out.kind = Proof::Timeout;
            return out;
        }
        int64_t target = A * x, y = ((target % C) * inv) % C;
        if (B * y > target) continue;
        int64_t z = (target - B * y) / C;
        if (A * x - B * y - C * z != 0) return out;
        out.kind = Proof::Optimal;
        out.x.assign(3, 0);
        out.x[v] = double(x);
        out.x[vars[0]] = double(y);
        out.x[vars[1]] = double(z);
        std::ostringstream js;
        js << "{\"type\":\"three_integer_congruence\",\"objective_col\":" << v << ",\"other_cols\":[" << vars[0] << ',' << vars[1]
           << "],\"inverse\":" << inv << ",\"start\":" << start << ",\"end\":" << x << '}';
        out.json = js.str();
        return out;
    }
    return out;
}
struct Grid {
    double step = 0, offset = 0;
    std::string json;
};
// Objective-only epigraph: minimize negative cost*t; t <= y_i and
// q*y_i = integral affine combination of integer columns. At an optimum,
// t equals min(y_i, finite upper bound), hence lies on 1/q grid. Merely
// having an integral-looking continuous cost is NOT enough.
inline Grid objective_grid(const Model& m) {
    Grid g;
    int n = int(m.cols.size()), nr = int(m.row_lo.size());
    if (n > 2048 || nr > 4096 || !m.qobj.empty() || !std::isfinite(m.obj_const)) return g;
    int t = -1;
    for (int j = 0; j < n; ++j)
        if (m.cost[j] != 0) {
            if (t >= 0) return g;
            t = j;
        }
    if (t < 0 || m.is_int[t] || !(m.cost[t] < 0) || m.col_lo[t] != 0) return g;
    std::vector<std::vector<Entry>> rows(nr);
    for (int j = 0; j < n; ++j)
        for (auto e : m.cols[j]) rows[e.index].push_back({j, e.value});
    if (m.cols[t].empty()) return g;
    int q = 0;
    std::ostringstream js;
    js << "{\"type\":\"objective_lattice\",\"objective_col\":" << t << ",\"links\":[";
    bool comma = false;
    for (auto e : m.cols[t]) {
        int i = e.index;
        if (m.row_lo[i] != -kInf || m.row_up[i] != 0 || e.value != 1 || rows[i].size() != 2) return g;
        int y = -1;
        for (auto a : rows[i])
            if (a.index != t) {
                if (a.value != -1) return g;
                y = a.index;
            }
        if (y < 0 || m.is_int[y] || m.col_lo[y] != 0 || std::isfinite(m.col_up[y])) return g;
        int er = -1, den = 0;
        for (auto a : m.cols[y]) {
            int r = a.index;
            if (r == i) continue;
            if (m.row_lo[r] != m.row_up[r] || !exact_int(m.row_lo[r]) || !exact_int(a.value) || a.value <= 0 || a.value > 1048576) return g;
            if (er >= 0) return g;
            er = r;
            den = int(a.value);
            for (auto v : rows[r])
                if (v.index != y && (!m.is_int[v.index] || !exact_int(v.value))) return g;
        }
        if (er < 0 || (q && den != q)) return g;
        q = den;
        if (comma) js << ',';
        comma = true;
        js << "{\"row\":" << i << ",\"grid_col\":" << y << ",\"equality\":" << er << '}';
    }
    if (std::isfinite(m.col_up[t]) && (!exact_int(m.col_up[t] * q) || m.col_up[t] < 0)) return g;
    g.step = -m.cost[t] / q;
    g.offset = m.obj_const;
    if (!std::isfinite(g.step) || g.step <= 0) return Grid{};
    js << "],\"denominator\":" << q << "}";
    g.json = js.str();
    return g;
}
inline double lattice_bound(double b, const Grid& g) {
    if (g.step == 0 || !std::isfinite(b)) return b;
    double units = (b - g.offset) / g.step;
    if (!std::isfinite(units) || std::abs(units) > 1099511627776.0) return b;
    // Conservative outward margin, never nearest rounding.
    double v = g.offset + std::ceil(units - 1e-7 * (1 + std::abs(units))) * g.step;
    return std::max(b, v);
}
} // namespace structural_integer
