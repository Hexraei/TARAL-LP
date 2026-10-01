// Sparse LU factorization written from first principles (no external library).
// Right-looking elimination on an active row-wise submatrix; pivots chosen by Markowitz
// cost (r-1)(c-1) among the sparsest few rows and columns, subject to a column threshold
// |a_pq| >= kThreshold * max_i |a_iq| that bounds the multipliers.
#include <algorithm>
#include <cmath>
#include <utility>

#include "taral.hpp"

namespace {
constexpr double kThreshold = 0.1;
constexpr double kTiny = 1e-11;  // below this a pivot is treated as zero (rank deficiency)
constexpr int kCandidates = 4;

double value_at(const std::vector<Entry>& row, int col) {
    for (const Entry& e : row)
        if (e.index == col) return e.value;
    return 0;
}
}  // namespace

bool SparseLU::factor(int m, const std::vector<std::vector<Entry>>& cols) {
    m_ = m;
    prow_.clear(), pcol_.clear(), piv_.clear(), L_.clear(), U_.clear();
    bad_pos.clear(), bad_rows.clear();

    std::vector<std::vector<Entry>> rows(m);  // active rows: (position, value)
    std::vector<std::vector<int>> colrows(m);
    for (int p = 0; p < m; ++p)
        for (const Entry& e : cols[p])
            if (e.value != 0) {
                rows[e.index].push_back({p, e.value});
                colrows[p].push_back(e.index);
            }
    std::vector<char> row_done(m, 0), col_done(m, 0);
    std::vector<int> where(m, -1);

    auto erase_from = [](std::vector<int>& v, int x) {
        for (size_t t = 0; t < v.size(); ++t)
            if (v[t] == x) {
                v[t] = v.back();
                v.pop_back();
                return;
            }
    };
    auto col_max = [&](int q) {
        double mx = 0;
        for (int i : colrows[q]) mx = std::max(mx, std::abs(value_at(rows[i], q)));
        return mx;
    };

    for (int k = 0; k < m; ++k) {
        // Sparsest few active columns and rows.
        int bc[kCandidates], br[kCandidates], nc = 0, nr = 0;
        auto keep = [](int* best, int& cnt, int id, size_t key, auto&& size_of) {
            if (cnt < kCandidates) best[cnt++] = id;
            else if (key < size_of(best[cnt - 1])) best[cnt - 1] = id;
            else return;
            for (int t = cnt - 1; t > 0 && size_of(best[t]) < size_of(best[t - 1]); --t) std::swap(best[t], best[t - 1]);
        };
        auto csize = [&](int q) { return colrows[q].size(); };
        auto rsize = [&](int i) { return rows[i].size(); };
        for (int q = 0; q < m; ++q)
            if (!col_done[q] && !colrows[q].empty()) keep(bc, nc, q, colrows[q].size(), csize);
        for (int i = 0; i < m; ++i)
            if (!row_done[i] && !rows[i].empty()) keep(br, nr, i, rows[i].size(), rsize);

        int p = -1, q = -1;
        double pv = 0, best_cost = 0;
        auto consider = [&](int i, int j, double v, double cmax) {
            if (std::abs(v) < kTiny || std::abs(v) < kThreshold * cmax) return;
            double cost = double(rows[i].size() - 1) * double(colrows[j].size() - 1);
            if (p < 0 || cost < best_cost || (cost == best_cost && std::abs(v) > std::abs(pv)))
                p = i, q = j, pv = v, best_cost = cost;
        };
        for (int t = 0; t < nc; ++t) {
            int j = bc[t];
            double cmax = col_max(j);
            for (int i : colrows[j]) consider(i, j, value_at(rows[i], j), cmax);
        }
        for (int t = 0; t < nr; ++t) {
            int i = br[t];
            for (const Entry& e : rows[i]) consider(i, e.index, e.value, col_max(e.index));
        }
        if (p < 0) {  // full scan before declaring rank deficiency
            for (int j = 0; j < m; ++j) {
                if (col_done[j]) continue;
                double cmax = col_max(j);
                for (int i : colrows[j]) consider(i, j, value_at(rows[i], j), cmax);
            }
        }
        if (p < 0) break;

        prow_.push_back(p), pcol_.push_back(q), piv_.push_back(pv);
        std::vector<Entry> urow;
        for (const Entry& e : rows[p])
            if (e.index != q) urow.push_back(e);
        U_.push_back(urow);

        std::vector<Entry> lcol;
        std::vector<int> targets = colrows[q];
        for (int i : targets) {
            if (i == p) continue;
            std::vector<Entry>& r = rows[i];
            for (size_t t = 0; t < r.size(); ++t) where[r[t].index] = int(t);
            double l = r[where[q]].value / pv;
            lcol.push_back({i, l});
            for (const Entry& e : urow) {
                if (where[e.index] >= 0) {
                    r[where[e.index]].value -= l * e.value;
                } else {
                    where[e.index] = int(r.size());
                    r.push_back({e.index, -l * e.value});
                    colrows[e.index].push_back(i);
                }
            }
            int at = where[q];
            for (const Entry& e : r) where[e.index] = -1;
            r[at] = r.back();
            r.pop_back();
        }
        L_.push_back(lcol);
        for (const Entry& e : urow) erase_from(colrows[e.index], p);
        colrows[q].clear();
        rows[p].clear();
        row_done[p] = col_done[q] = 1;
    }

    if (int(prow_.size()) == m) return true;
    for (int j = 0; j < m; ++j)
        if (!col_done[j]) bad_pos.push_back(j);
    for (int i = 0; i < m; ++i)
        if (!row_done[i]) bad_rows.push_back(i);
    return false;
}

void SparseLU::ftran(std::vector<double>& rhs, std::vector<double>& out) const {
    for (size_t k = 0; k < L_.size(); ++k) {
        double v = rhs[prow_[k]];
        if (v != 0)
            for (const Entry& e : L_[k]) rhs[e.index] -= e.value * v;
    }
    out.assign(m_, 0);
    for (int k = m_ - 1; k >= 0; --k) {
        double s = rhs[prow_[k]];
        for (const Entry& e : U_[k]) s -= e.value * out[e.index];
        out[pcol_[k]] = s / piv_[k];
    }
}

void SparseLU::btran(std::vector<double>& rhs, std::vector<double>& out) const {
    out.assign(m_, 0);
    for (int k = 0; k < m_; ++k) {
        double w = rhs[pcol_[k]] / piv_[k];
        out[prow_[k]] = w;
        if (w != 0)
            for (const Entry& e : U_[k]) rhs[e.index] -= e.value * w;
    }
    for (int k = m_ - 1; k >= 0; --k) {
        double s = 0;
        for (const Entry& e : L_[k]) s += e.value * out[e.index];
        out[prow_[k]] -= s;
    }
}
