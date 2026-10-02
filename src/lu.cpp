// Sparse LU factorization written from first principles (no external library).
// Right-looking elimination on the active submatrix. Singleton columns and rows pivot
// directly; the remaining nucleus uses Markowitz cost (r-1)(c-1) among the sparsest few
// rows and columns, subject to a column threshold |a_pq| >= kThreshold * max_i |a_iq|
// that bounds the multipliers.
#include <algorithm>
#include <cmath>
#include <set>
#include <utility>

#include "taral.hpp"

namespace {
constexpr double kThreshold = 0.1;
constexpr double kTiny = 1e-11;  // below this a pivot is treated as zero (rank deficiency)
constexpr int kCandidates = 4;

// One matrix entry as stored in the active submatrix. It exists twice, once in its row (index = column)
// and once in its column (index = row); `pos` is where the twin sits in the other list, so an update or a
// removal never searches a list.
struct Cell {
    int index;
    double value;
    int pos;
};
}  // namespace

bool SparseLU::factor(int m, const std::vector<std::vector<Entry>>& cols) {
    m_ = m;
    prow_.clear(), pcol_.clear(), piv_.clear(), L_.clear(), U_.clear();
    bad_pos.clear(), bad_rows.clear();

    // The active submatrix is kept both row-wise (position, value) and column-wise
    // (row, value), so pivot search reads column values without scanning long rows.
    std::vector<std::vector<Cell>> rows(m), acol(m);
    for (int p = 0; p < m; ++p)
        for (const Entry& e : cols[p])
            if (e.value != 0) {
                int rpos = int(rows[e.index].size()), cpos = int(acol[p].size());
                rows[e.index].push_back({p, e.value, cpos});
                acol[p].push_back({e.index, e.value, rpos});
            }
    std::vector<char> row_done(m, 0), col_done(m, 0);
    std::vector<int> where(m, -1);
    std::vector<int> col_single, row_single;  // candidates that may have exactly one entry
    for (int j = 0; j < m; ++j)
        if (acol[j].size() == 1) col_single.push_back(j);
    for (int i = 0; i < m; ++i)
        if (rows[i].size() == 1) row_single.push_back(i);
    // Active nonempty columns and rows ordered by (count, index): the Markowitz search takes the
    // first few, the same choice (ties to the lower index) as a full scan, without the O(m) pass.
    std::set<std::pair<size_t, int>> cset, rset;
    for (int j = 0; j < m; ++j)
        if (!acol[j].empty()) cset.insert({acol[j].size(), j});
    for (int i = 0; i < m; ++i)
        if (!rows[i].empty()) rset.insert({rows[i].size(), i});

    std::vector<double> cmax(m, -1);  // cached max |value| per active column, -1 = stale
    auto col_max = [&](int q) {
        if (cmax[q] < 0) {
            double mx = 0;
            for (const Cell& e : acol[q]) mx = std::max(mx, std::abs(e.value));
            cmax[q] = mx;
        }
        return cmax[q];
    };

    for (int k = 0; k < m; ++k) {
        int p = -1, q = -1;
        double pv = 0, best_cost = 0;
        auto consider = [&](int i, int j, double v, double cmax) {
            if (std::abs(v) < kTiny || std::abs(v) < kThreshold * cmax) return;
            double cost = double(rows[i].size() - 1) * double(acol[j].size() - 1);
            if (p < 0 || cost < best_cost || (cost == best_cost && std::abs(v) > std::abs(pv)))
                p = i, q = j, pv = v, best_cost = cost;
        };
        auto consider_col = [&](int j) {
            double cmax = col_max(j);
            for (const Cell& e : acol[j]) consider(e.index, j, e.value, cmax);
        };

        while (p < 0 && !col_single.empty()) {
            int j = col_single.back();
            col_single.pop_back();
            if (!col_done[j] && acol[j].size() == 1) consider_col(j);
        }
        while (p < 0 && !row_single.empty()) {
            int i = row_single.back();
            row_single.pop_back();
            if (!row_done[i] && rows[i].size() == 1) consider(i, rows[i][0].index, rows[i][0].value, col_max(rows[i][0].index));
        }
        if (p < 0) {  // Markowitz search among the sparsest few active columns and rows
            int bc[kCandidates], br[kCandidates], nc = 0, nr = 0;
            for (auto it = cset.begin(); it != cset.end() && nc < kCandidates; ++it) bc[nc++] = it->second;
            for (auto it = rset.begin(); it != rset.end() && nr < kCandidates; ++it) br[nr++] = it->second;
            for (int t = 0; t < nc; ++t) consider_col(bc[t]);
            for (int t = 0; t < nr; ++t)
                for (const Cell& e : rows[br[t]]) consider(br[t], e.index, e.value, col_max(e.index));
            if (p < 0)  // full scan before declaring rank deficiency
                for (int j = 0; j < m; ++j)
                    if (!col_done[j]) consider_col(j);
        }
        if (p < 0) break;

        prow_.push_back(p), pcol_.push_back(q), piv_.push_back(pv);
        // Only the pivot row/column and the rows of column q / columns of row p change counts.
        std::vector<int> mrows, mcols;
        for (const Cell& t : acol[q]) mrows.push_back(t.index), rset.erase({rows[t.index].size(), t.index});
        for (const Cell& e : rows[p]) mcols.push_back(e.index), cset.erase({acol[e.index].size(), e.index});
        std::vector<Cell> urow;  // pos = where row p sits in that column's list
        for (const Cell& e : rows[p])
            if (e.index != q) urow.push_back(e);
        U_.emplace_back();
        U_.back().reserve(urow.size());
        for (const Cell& e : urow) U_.back().push_back({e.index, e.value});

        std::vector<Entry> lcol;
        lcol.reserve(acol[q].size());
        for (const Cell& t : acol[q]) {
            int i = t.index;
            if (i == p) continue;
            std::vector<Cell>& r = rows[i];
            for (size_t s = 0; s < r.size(); ++s) where[r[s].index] = int(s);
            double l = t.value / pv;
            lcol.push_back({i, l});
            for (const Cell& e : urow) {
                double delta = -l * e.value;
                if (where[e.index] >= 0) {
                    Cell& c = r[where[e.index]];
                    c.value += delta;
                    acol[e.index][c.pos].value += delta;
                } else {
                    where[e.index] = int(r.size());
                    std::vector<Cell>& col = acol[e.index];
                    r.push_back({e.index, delta, int(col.size())});
                    col.push_back({i, delta, int(r.size()) - 1});
                }
            }
            int at = where[q];
            for (const Cell& e : r) where[e.index] = -1;
            r[at] = r.back();
            r.pop_back();
            if (at < int(r.size())) acol[r[at].index][r[at].pos].pos = at;  // the moved entry's twin
            if (r.size() == 1) row_single.push_back(i);
        }
        L_.push_back(lcol);
        for (const Cell& e : urow) {
            cmax[e.index] = -1;  // values changed by the update above, and row p leaves the column
            std::vector<Cell>& col = acol[e.index];
            col[e.pos] = col.back();
            col.pop_back();
            if (e.pos < int(col.size())) rows[col[e.pos].index][col[e.pos].pos].pos = e.pos;
            if (col.size() == 1) col_single.push_back(e.index);
        }
        acol[q].clear();
        rows[p].clear();
        row_done[p] = col_done[q] = 1;
        for (int i : mrows)
            if (!row_done[i] && !rows[i].empty()) rset.insert({rows[i].size(), i});
        for (int j : mcols)
            if (!col_done[j] && !acol[j].empty()) cset.insert({acol[j].size(), j});
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
