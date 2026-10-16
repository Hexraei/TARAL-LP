#pragma GCC optimize("O3")
#pragma GCC optimize("fp-contract=off")
#include "ipm_super.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
namespace taral_qp_super {
using std::vector; using std::max; using std::min; using std::fabs;
struct SuperImpl {
    int n = 0, ns = 0; vector<int> Kp, Ki; vector<signed char> sg;
    vector<int> sbeg, sparent, rowp, rows; // supernode s: cols [sbeg[s], sbeg[s+1]); rows (sorted, first width are its cols) in rows[rowp[s]..rowp[s+1])
    vector<long long> lp; // offset of dense panel (f x w col-major, ld=f)
    Vec Lx, D; vector<vector<int>> children; long nreg = 0; double maxfront = 0; long long fl = 0;
    // column-wise lower structure of A (CSC): for col j rows >= j
    vector<long long> Ap; vector<int> Ai; vector<long long> Amap;
    int relax = 0;
    void analyze(int n_, const vector<int>& Kp_, const vector<int>& Ki_, const vector<signed char>& sg_, const vector<int>& parent, const vector<int>& cnt) {
        n = n_; Kp = Kp_; Ki = Ki_; sg = sg_;
        // supernode partition
        vector<int> start{0};
        for (int j = 1; j < n; ++j) { bool merge = parent[j - 1] == j && (cnt[j - 1] - 1 - cnt[j]) <= relax && (cnt[j - 1] - 1 - cnt[j]) >= 0; if (!merge) start.push_back(j); }
        start.push_back(n); ns = int(start.size()) - 1; sbeg = start;
        vector<int> snode(n); for (int s = 0; s < ns; ++s) for (int j = sbeg[s]; j < sbeg[s + 1]; ++j) snode[j] = s;
        sparent.assign(ns, -1); children.assign(ns, {});
        for (int s = 0; s < ns; ++s) { int last = sbeg[s + 1] - 1; if (parent[last] >= 0) { sparent[s] = snode[parent[last]]; children[sparent[s]].push_back(s); } }
        // CSC lower of A
        vector<long long> c(n + 1, 0); for (int k = 0; k < n; ++k) for (int p = Kp[k]; p < Kp[k + 1]; ++p) ++c[Ki[p] + 1];
        Ap.assign(n + 1, 0); for (int j = 0; j < n; ++j) Ap[j + 1] = Ap[j] + c[j + 1]; Ai.resize(Ap[n]); Amap.resize(Ap[n]); vector<long long> pos(Ap.begin(), Ap.end() - 1);
        for (int k = 0; k < n; ++k) for (int p = Kp[k]; p < Kp[k + 1]; ++p) { long long q = pos[Ki[p]]++; Ai[q] = k; Amap[q] = p; }
        // row structures
        rowp.assign(ns + 1, 0); rows.clear(); vector<int> mark(n, -1);
        for (int s = 0; s < ns; ++s) {
            int c0 = sbeg[s], c1 = sbeg[s + 1]; rowp[s] = int(rows.size()); size_t b0 = rows.size();
            for (int j = c0; j < c1; ++j) { mark[j] = s; rows.push_back(j); }
            for (int j = c0; j < c1; ++j) for (long long q = Ap[j]; q < Ap[j + 1]; ++q) { int r = Ai[q]; if (r >= c1 && mark[r] != s) { mark[r] = s; rows.push_back(r); } }
            for (int ch : children[s]) for (int q = rowp[ch]; q < rowp[ch + 1]; ++q) { int r = rows[q]; if (r >= c1 && mark[r] != s) { mark[r] = s; rows.push_back(r); } }
            sort(rows.begin() + b0 + (c1 - c0), rows.end()); rowp[s + 1] = int(rows.size());
            // rowp for later reads of child uses rowp[ch+1]: set before use (children precede parents)
        }
        lp.assign(ns + 1, 0); fl = 0;
        for (int s = 0; s < ns; ++s) { long long w = sbeg[s + 1] - sbeg[s], f = rowp[s + 1] - rowp[s]; lp[s + 1] = lp[s] + w * f; maxfront = max(maxfront, double(f)); for (long long k = 0; k < w; ++k) fl += (f - k) * (f - k); }
        Lx.assign(lp[ns], 0); D.assign(n, 0);
    }
    // returns false when the deadline passes
    __attribute__((target_clones("avx2", "default")))
    bool factor(const Vec& V, Clock::time_point deadline, double kPivTol, double kPivRep, double kPivAbs) {
        nreg = 0; vector<Vec> U(ns), BU(ns); vector<int> relpos(n, -1);
        Vec F;
        for (int s = 0; s < ns; ++s) {
            if ((s & 63) == 0 && Clock::now() > deadline) return false;
            int c0 = sbeg[s], w = sbeg[s + 1] - c0, f = rowp[s + 1] - rowp[s]; const int* R = &rows[rowp[s]];
            for (int i = 0; i < f; ++i) relpos[R[i]] = i;
            F.assign((size_t)f * f, 0.0); // use lower triangle, col-major ld=f
            Vec dg(w, 0.0); Vec bg(f, 0.0);
            for (int j = c0; j < sbeg[s + 1]; ++j) for (long long q = Ap[j]; q < Ap[j + 1]; ++q) { int r = relpos[Ai[q]]; double v = V[Amap[q]]; F[(size_t)(j - c0) * f + r] += v; if (r == j - c0) dg[j - c0] = fabs(v); }
            for (int j = 0; j < w; ++j) bg[j] = dg[j];
            for (int ch : children[s]) {
                int cw = sbeg[ch + 1] - sbeg[ch], cf = rowp[ch + 1] - rowp[ch], nb = cf - cw; const int* RC = &rows[rowp[ch]] + cw; Vec& u = U[ch];
                for (int b = 0; b < nb; ++b) { int cb = relpos[RC[b]]; for (int a = b; a < nb; ++a) { int ra = relpos[RC[a]]; F[(size_t)cb * f + ra] += u[(size_t)b * nb + a]; } }
                Vec().swap(u);
                { Vec& bu = BU[ch]; for (int b = 0; b < nb; ++b) { int cb = relpos[RC[b]]; bg[cb] = max(bg[cb], bu[b]); } Vec().swap(bu); }
            }
            // dense partial LDL^T of the first w columns, panel width PB, fused update of 8 columns
            double* Fp = F.data(); const int PB = 32;
            Vec dd(w);
            for (int kb = 0; kb < w; kb += PB) {
                int kend = min(w, kb + PB);
                for (int k = kb; k < kend; ++k) {
                    // update column k with panel columns kb..k-1 (rows k..f)
                    for (int m = kb; m < k; ++m) { double t = dd[m] * Fp[(size_t)m * f + k]; const double* Lm = Fp + (size_t)m * f; double* Ck = Fp + (size_t)k * f; for (int i = k; i < f; ++i) Ck[i] -= Lm[i] * t; }
                    double d = Fp[(size_t)k * f + k]; double big = max(bg[k], fabs(d));
                    int gk = c0 + k;
                    if (!(sg[gk] * d > kPivTol * big)) d = sg[gk] * max(kPivAbs, kPivRep * big), ++nreg;
                    dd[k] = d; double inv = 1.0 / d; double* Ck = Fp + (size_t)k * f; for (int i = k + 1; i < f; ++i) Ck[i] *= inv; Ck[k] = d; for (int i = k + 1; i < f; ++i) bg[i] = max(bg[i], fabs(Ck[i] * Ck[i] * d)); // store L (unit diag implicit), D on diagonal
                }
                // trailing update: columns j in [kend, f) get -= sum_{m in panel} L[:,m] * dd[m] * L[j,m]
                for (int j = kend; j < f; ++j) {
                    double* Cj = Fp + (size_t)j * f; int m = kb;
                    for (; m + 8 <= kend; m += 8) {
                        double t0 = dd[m] * Fp[(size_t)m * f + j], t1 = dd[m + 1] * Fp[(size_t)(m + 1) * f + j], t2 = dd[m + 2] * Fp[(size_t)(m + 2) * f + j], t3 = dd[m + 3] * Fp[(size_t)(m + 3) * f + j];
                        double t4 = dd[m + 4] * Fp[(size_t)(m + 4) * f + j], t5 = dd[m + 5] * Fp[(size_t)(m + 5) * f + j], t6 = dd[m + 6] * Fp[(size_t)(m + 6) * f + j], t7 = dd[m + 7] * Fp[(size_t)(m + 7) * f + j];
                        const double *L0 = Fp + (size_t)m * f, *L1 = L0 + f, *L2 = L1 + f, *L3 = L2 + f, *L4 = L3 + f, *L5 = L4 + f, *L6 = L5 + f, *L7 = L6 + f;
                        for (int i = j; i < f; ++i) Cj[i] -= L0[i] * t0 + L1[i] * t1 + L2[i] * t2 + L3[i] * t3 + L4[i] * t4 + L5[i] * t5 + L6[i] * t6 + L7[i] * t7;
                    }
                    for (; m < kend; ++m) { double t = dd[m] * Fp[(size_t)m * f + j]; const double* Lm = Fp + (size_t)m * f; for (int i = j; i < f; ++i) Cj[i] -= Lm[i] * t; }
                }
            }
            // store panel and D
            double* out = &Lx[lp[s]]; for (int k = 0; k < w; ++k) { memcpy(out + (size_t)k * f, Fp + (size_t)k * f, sizeof(double) * f); D[c0 + k] = dd[k]; }
            int nb = f - w;
            if (nb > 0 && sparent[s] >= 0) { Vec& u = U[s]; u.assign((size_t)nb * nb, 0.0); for (int b = 0; b < nb; ++b) memcpy(&u[(size_t)b * nb + b], Fp + (size_t)(w + b) * f + (w + b), sizeof(double) * (nb - b));
                Vec& bu = BU[s]; bu.assign(bg.begin() + w, bg.end()); }
        }
        return true;
    }
    void solve(Vec& x) const {
        Vec t;
        for (int s = 0; s < ns; ++s) { int c0 = sbeg[s], w = sbeg[s + 1] - c0, f = rowp[s + 1] - rowp[s]; const int* R = &rows[rowp[s]]; const double* L = &Lx[lp[s]];
            for (int k = 0; k < w; ++k) { double xk = x[c0 + k]; const double* Lk = L + (size_t)k * f; for (int i = k + 1; i < f; ++i) x[R[i]] -= Lk[i] * xk; } }
        for (int j = 0; j < n; ++j) x[j] /= D[j];
        for (int s = ns - 1; s >= 0; --s) { int c0 = sbeg[s], w = sbeg[s + 1] - c0, f = rowp[s + 1] - rowp[s]; const int* R = &rows[rowp[s]]; const double* L = &Lx[lp[s]];
            for (int k = w - 1; k >= 0; --k) { double sacc = x[c0 + k]; const double* Lk = L + (size_t)k * f; for (int i = k + 1; i < f; ++i) sacc -= Lk[i] * x[R[i]]; x[c0 + k] = sacc; } }
    }
};

SuperLdl::SuperLdl() : impl_(new SuperImpl) {}
SuperLdl::~SuperLdl() { delete impl_; }
void SuperLdl::analyze(int n, const vector<int>& Kp, const vector<int>& Ki, const vector<signed char>& sign, const vector<int>& parent, const vector<int>& cnt) { impl_->analyze(n, Kp, Ki, sign, parent, cnt); }
bool SuperLdl::factor(const Vec& V, Clock::time_point deadline, double a, double b, double c) { return impl_->factor(V, deadline, a, b, c); }
void SuperLdl::solve(Vec& x) const { impl_->solve(x); }
bool SuperLdl::all_finite() const { for (double v : impl_->D) if (!std::isfinite(v)) return false; for (double v : impl_->Lx) if (!std::isfinite(v)) return false; return true; }
int SuperLdl::nreg() const { return int(impl_->nreg); }
vector<int> etree_postorder(const vector<int>& parent) {
    int n = int(parent.size()); vector<int> head(n, -1), nxt(n, -1);
    for (int j = n - 1; j >= 0; --j) if (parent[j] >= 0) { nxt[j] = head[parent[j]]; head[parent[j]] = j; }
    vector<int> post; post.reserve(n); vector<int> stk;
    for (int r = 0; r < n; ++r) if (parent[r] < 0) { stk.push_back(r); while (!stk.empty()) { int v = stk.back(); if (head[v] >= 0) { int c = head[v]; head[v] = nxt[c]; stk.push_back(c); } else { post.push_back(v); stk.pop_back(); } } }
    return post;
}
}
