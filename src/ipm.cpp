// Primal-dual interior-point method (Mehrotra predictor-corrector) for LP and convex QP:
//   min c'x + 0.5 x'Qx  s.t.  row_lo <= Ax <= row_up,  col_lo <= x <= col_up.
// Pipeline: light presolve (fixed columns substituted, empty/free rows dropped, empty columns set
// directly), convexity check of Q by an LDL^T inertia test, geometric + Ruiz scaling, then the IPM on
//   min c'z + 0.5 z'Qz  s.t.  Az = b,  l <= z <= u   (one slack per inequality / ranged row).
// Each finite bound gets its own slack/dual pair (z - xl = l, z + xu = u); free variables have none.
// Newton systems are solved on the quasi-definite augmented matrix
//   [ -(Q + Theta^-1 + rho I)  A' ;  A  delta I ]
// with an own approximate-minimum-degree ordering, an own up-looking sparse LDL^T with static and
// dynamic regularisation, and iterative refinement against the unregularised matrix.
// Infeasibility / unboundedness: no homogeneous embedding. Every iteration the current dual iterate y
// and the last dual direction are tested as Farkas certificates, and the primal iterate and last primal
// direction as recession directions (A d = 0, Q d = 0, c'd < 0, d inside the bound cone). A status is
// only reported when such a certificate passes the tolerances below.
// "Optimal" requires the original-model measures (recomputed from x and y alone) below tol.
#include "ipm.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>

namespace {

using Vec = std::vector<double>;
using Clock = std::chrono::steady_clock;
constexpr double kBig = 1e20;
const double kRho = std::getenv("IPM_RHO") ? std::atof(std::getenv("IPM_RHO")) : 1e-8;
const double kDelta = std::getenv("IPM_DELTA") ? std::atof(std::getenv("IPM_DELTA")) : 1e-8;  // |bound| >= kBig counts as infinite (the HiGHS convention)
const double kSlow = std::getenv("IPM_SLOW") ? std::atof(std::getenv("IPM_SLOW")) : 0.9;
const double kPivTol = std::getenv("IPM_PIV") ? std::atof(std::getenv("IPM_PIV")) : 1e-30;
const double kPivRep = std::getenv("IPM_PIVREP") ? std::atof(std::getenv("IPM_PIVREP")) : 1e-10;
const double kPivAbs = std::getenv("IPM_PIVABS") ? std::atof(std::getenv("IPM_PIVABS")) : 2e-7;
const double kKeepDelta = std::getenv("IPM_KEEPD") ? std::atof(std::getenv("IPM_KEEPD")) : 1;
const int kRefMax = std::getenv("IPM_REF") ? std::atoi(std::getenv("IPM_REF")) : 30;

double norm_inf(const Vec& v) {
    double m = 0;
    for (double a : v) m = std::max(m, std::abs(a));
    return m;
}

// ---------------------------------------------------------------------------------------------
// Approximate minimum degree ordering on a quotient graph: element absorption (including aggressive
// absorption) and the approximate external degree bound of Amestoy, Davis and Duff. No supervariable
// detection or mass elimination. Rows/columns with more than 10 sqrt(n) neighbours are ordered last.
// adj must be symmetric without self loops. Returns order[k] = original index of the k-th pivot.
// ponytail: no supervariables; add them if ordering time shows up on large problems.
std::vector<int> amd_order(int n, std::vector<std::vector<int>> adj) {
    std::vector<std::vector<int>> E(n), L(n);
    std::vector<char> st(n, 0);  // 0 variable, 1 element, 2 absorbed element, 3 dense
    const size_t dense = std::max<size_t>(16, static_cast<size_t>(10 * std::sqrt(static_cast<double>(n))));
    std::vector<int> order, dense_nodes;
    order.reserve(n);
    for (int i = 0; i < n; ++i)
        if (adj[i].size() > dense) st[i] = 3, dense_nodes.push_back(i);
    std::vector<int> deg(n, 0), head(n + 1, -1), nxt(n, -1), prv(n, -1);
    auto ins = [&](int i) {
        int d = deg[i];
        nxt[i] = head[d], prv[i] = -1;
        if (head[d] >= 0) prv[head[d]] = i;
        head[d] = i;
    };
    auto del = [&](int i) {
        if (prv[i] >= 0) nxt[prv[i]] = nxt[i];
        else head[deg[i]] = nxt[i];
        if (nxt[i] >= 0) prv[nxt[i]] = prv[i];
    };
    int nleft = 0;
    for (int i = 0; i < n; ++i) {
        if (st[i]) continue;
        auto& a = adj[i];
        a.erase(std::remove_if(a.begin(), a.end(), [&](int v) { return st[v] == 3; }), a.end());
        deg[i] = static_cast<int>(a.size());
        ins(i);
        ++nleft;
    }
    std::vector<int> mark(n, 0), wflag(n, 0), w(n, 0), Lp;
    int stamp = 0, wstamp = 0, mindeg = 0;
    while (nleft > 0) {
        while (head[mindeg] < 0) ++mindeg;
        const int p = head[mindeg];
        del(p);
        mark[p] = ++stamp;
        Lp.clear();
        for (int e : E[p]) {
            if (st[e] != 1) continue;
            for (int v : L[e])
                if (st[v] == 0 && mark[v] != stamp) mark[v] = stamp, Lp.push_back(v);
            st[e] = 2;
            std::vector<int>().swap(L[e]);
        }
        for (int v : adj[p])
            if (st[v] == 0 && mark[v] != stamp) mark[v] = stamp, Lp.push_back(v);
        std::vector<int>().swap(adj[p]);
        std::vector<int>().swap(E[p]);
        st[p] = 1;
        order.push_back(p);
        --nleft;
        ++wstamp;  // w[e] = |L_e \ L_p| for elements adjacent to L_p
        for (int i : Lp)
            for (int e : E[i]) {
                if (st[e] != 1) continue;
                if (wflag[e] != wstamp) wflag[e] = wstamp, w[e] = static_cast<int>(L[e].size());
                --w[e];
            }
        const long lp = static_cast<long>(Lp.size());
        for (int i : Lp) {
            del(i);
            long ext = 0;
            size_t k = 0;
            for (int e : E[i]) {
                if (st[e] != 1) continue;
                if (w[e] == 0) {  // L_e inside L_p: absorb e into p
                    st[e] = 2;
                    std::vector<int>().swap(L[e]);
                    continue;
                }
                E[i][k++] = e;
                ext += w[e];
            }
            E[i].resize(k);
            E[i].push_back(p);
            k = 0;
            for (int v : adj[i])
                if (st[v] == 0 && mark[v] != stamp) adj[i][k++] = v;
            adj[i].resize(k);
            long d = static_cast<long>(adj[i].size()) + (lp - 1) + ext;
            d = std::min({d, static_cast<long>(deg[i]) + lp - 1, static_cast<long>(nleft) - 1});
            deg[i] = static_cast<int>(std::max(d, 0L));
            ins(i);
            mindeg = std::min(mindeg, deg[i]);
        }
        L[p] = Lp;
    }
    std::sort(dense_nodes.begin(), dense_nodes.end(),
              [&](int a, int b) { return adj[a].size() < adj[b].size() || (adj[a].size() == adj[b].size() && a < b); });
    order.insert(order.end(), dense_nodes.begin(), dense_nodes.end());
    return order;
}

// ---------------------------------------------------------------------------------------------
// Sparse LDL^T (up-looking, elimination tree; after Davis' LDL) of a symmetric matrix given as a
// list of entries (one triangle plus the full diagonal). Values come from V[src[e]], so the same
// symbolic analysis serves every refactorisation. Each pivot has an expected sign; with
// regularisation on, a pivot whose sign is wrong or magnitude below dyn_eps is replaced by
// sign * dyn_delta (quasi-definite systems are strongly factorisable in any symmetric order).
class Ldl {
public:
    void analyze(int n, const std::vector<int>& ei, const std::vector<int>& ej, const std::vector<int>& esrc,
                 const std::vector<signed char>& sign) {
        n_ = n;
        std::vector<std::vector<int>> adj(n);
        for (size_t e = 0; e < ei.size(); ++e)
            if (ei[e] != ej[e]) adj[ei[e]].push_back(ej[e]), adj[ej[e]].push_back(ei[e]);
        perm_ = amd_order(n, std::move(adj));
        std::vector<int> iperm(n);
        for (int k = 0; k < n; ++k) iperm[perm_[k]] = k;
        Kp_.assign(n + 1, 0);
        for (size_t e = 0; e < ei.size(); ++e) ++Kp_[std::max(iperm[ei[e]], iperm[ej[e]]) + 1];
        for (int k = 0; k < n; ++k) Kp_[k + 1] += Kp_[k];
        Ki_.resize(ei.size());
        Ksrc_.resize(ei.size());
        std::vector<int> pos(Kp_.begin(), Kp_.end() - 1);
        for (size_t e = 0; e < ei.size(); ++e) {
            int a = iperm[ei[e]], b = iperm[ej[e]];
            int q = pos[std::max(a, b)]++;
            Ki_[q] = std::min(a, b);
            Ksrc_[q] = esrc[e];
        }
        sign_.resize(n);
        for (int k = 0; k < n; ++k) sign_[k] = sign[perm_[k]];
        parent_.assign(n, -1);
        flag_.assign(n, -1);
        lnz_.assign(n, 0);
        for (int k = 0; k < n; ++k) {
            flag_[k] = k;
            for (int p = Kp_[k]; p < Kp_[k + 1]; ++p)
                for (int i = Ki_[p]; flag_[i] != k; i = parent_[i]) {
                    if (parent_[i] == -1) parent_[i] = k;
                    ++lnz_[i];
                    flag_[i] = k;
                }
        }
        Lp_.assign(n + 1, 0);
        for (int k = 0; k < n; ++k) Lp_[k + 1] = Lp_[k] + lnz_[k];
        y_.assign(n, 0);
        pattern_.assign(n, 0);
        D_.assign(n, 0);
        x_.assign(n, 0);
    }

    long long nnz_l() const { return n_ ? Lp_[n_] : 0; }
    bool fits() const { return nnz_l() <= 400000000LL; }

    // Returns false when the deadline passes.
    bool factor(const Vec& V, Clock::time_point deadline, bool regularize) {
        if (Li_.size() != static_cast<size_t>(nnz_l())) Li_.resize(nnz_l()), Lx_.resize(nnz_l());
        nreg = 0, nwrong = 0;
        for (int k = 0; k < n_; ++k) {
            if ((k & 255) == 0 && Clock::now() > deadline) return false;
            y_[k] = 0;
            int top = n_;
            flag_[k] = k;
            lnz_[k] = 0;
            for (int p = Kp_[k]; p < Kp_[k + 1]; ++p) {
                int i = Ki_[p];
                y_[i] += V[Ksrc_[p]];
                int len = 0;
                for (; flag_[i] != k; i = parent_[i]) pattern_[len++] = i, flag_[i] = k;
                while (len > 0) pattern_[--top] = pattern_[--len];
            }
            double d = y_[k], big = std::abs(d);
            y_[k] = 0;
            for (; top < n_; ++top) {
                const int i = pattern_[top];
                const double yi = y_[i];
                y_[i] = 0;
                const long long p0 = Lp_[i], p2 = Lp_[i] + lnz_[i];
                for (long long p = p0; p < p2; ++p) y_[Li_[p]] -= Lx_[p] * yi;
                const double lki = yi / D_[i];
                d -= lki * yi;
                big = std::max(big, std::abs(lki * yi));
                Li_[p2] = k;
                Lx_[p2] = lki;
                ++lnz_[i];
            }
            if (sign_[k] * d <= 0) ++nwrong;
            // A pivot of the wrong sign or lost in cancellation (relative to the terms that formed it)
            // is replaced by a huge value, which decouples that row (the Krylov solve then corrects).
            if (regularize && !(sign_[k] * d > kPivTol * big)) d = sign_[k] * std::max(kPivAbs, kPivRep * big), ++nreg;
            D_[k] = d;
        }
        return true;
    }

    void solve(Vec& b) const {
        for (int k = 0; k < n_; ++k) x_[k] = b[perm_[k]];
        for (int k = 0; k < n_; ++k) {
            const double xk = x_[k];
            if (xk != 0)
                for (long long p = Lp_[k]; p < Lp_[k + 1]; ++p) x_[Li_[p]] -= Lx_[p] * xk;
        }
        for (int k = 0; k < n_; ++k) x_[k] /= D_[k];
        for (int k = n_ - 1; k >= 0; --k) {
            double s = x_[k];
            for (long long p = Lp_[k]; p < Lp_[k + 1]; ++p) s -= Lx_[p] * x_[Li_[p]];
            x_[k] = s;
        }
        for (int k = 0; k < n_; ++k) b[perm_[k]] = x_[k];
    }

    int nreg = 0, nwrong = 0;

private:
    int n_ = 0;
    std::vector<int> perm_, Kp_, Ki_, Ksrc_, parent_, flag_, lnz_, pattern_, Li_;
    std::vector<signed char> sign_;
    std::vector<long long> Lp_;
    Vec Lx_, D_, y_;
    mutable Vec x_;
};

// ---------------------------------------------------------------------------------------------
// Internal problem: min c'z + 0.5 z'Qz  s.t.  A z = b,  l <= z <= u  (infinite bounds are +-inf).
struct Qp {
    int N = 0, M = 0;
    std::vector<int> Ap, Ai;  // CSC, M x N
    Vec Ax;
    Vec qd;                   // diagonal of Q
    std::vector<int> qr, qc;  // strictly lower off-diagonal entries (qr > qc)
    Vec qv;
    Vec c, b, l, u;

    void mul_A(const Vec& z, Vec& out) const {
        out.assign(M, 0);
        for (int j = 0; j < N; ++j)
            if (z[j] != 0)
                for (int p = Ap[j]; p < Ap[j + 1]; ++p) out[Ai[p]] += Ax[p] * z[j];
    }
    void mul_At(const Vec& y, Vec& out) const {
        out.assign(N, 0);
        for (int j = 0; j < N; ++j) {
            double s = 0;
            for (int p = Ap[j]; p < Ap[j + 1]; ++p) s += Ax[p] * y[Ai[p]];
            out[j] = s;
        }
    }
    void mul_Q(const Vec& z, Vec& out) const {
        out.resize(N);
        for (int j = 0; j < N; ++j) out[j] = qd[j] * z[j];
        for (size_t k = 0; k < qv.size(); ++k) out[qr[k]] += qv[k] * z[qc[k]], out[qc[k]] += qv[k] * z[qr[k]];
    }
};

struct RunOut {
    Vec z, y;
    IpmStatus status = IpmStatus::NumericalFailure;
    long iters = 0;
    std::string msg;
    long long fnnz = 0;
    Vec ray;  // recession direction behind a DualInfeasible status (internal scale, max-norm 1)
};

// accept(z, y) re-checks a candidate on the original model and returns true when it passes.
RunOut run_ipm(const Qp& P, const IpmOptions& opt, Clock::time_point deadline,
               const std::function<bool(const Vec&, const Vec&)>& accept) {
    const int N = P.N, M = P.M;
    RunOut out;
    std::vector<char> hl(N), hu(N);
    int nc = 0;
    for (int j = 0; j < N; ++j) {
        hl[j] = std::isfinite(P.l[j]);
        hu[j] = std::isfinite(P.u[j]);
        nc += hl[j] + hu[j];
    }
    const int nqo = static_cast<int>(P.qv.size()), nza = P.Ap[N];
    // Augmented system pattern and value sources.
    std::vector<int> ei, ej, es;
    std::vector<signed char> sign(N + M);
    for (int j = 0; j < N; ++j) ei.push_back(j), ej.push_back(j), es.push_back(j), sign[j] = -1;
    for (int k = 0; k < nqo; ++k) ei.push_back(P.qr[k]), ej.push_back(P.qc[k]), es.push_back(N + k);
    for (int j = 0; j < N; ++j)
        for (int p = P.Ap[j]; p < P.Ap[j + 1]; ++p) ei.push_back(N + P.Ai[p]), ej.push_back(j), es.push_back(N + nqo + p);
    for (int i = 0; i < M; ++i) ei.push_back(N + i), ej.push_back(N + i), es.push_back(N + nqo + nza + i), sign[N + i] = 1;
    Ldl ldl;
    ldl.analyze(N + M, ei, ej, es, sign);
    out.fnnz = ldl.nnz_l();
    if (!ldl.fits()) {
        out.msg = "augmented-system factor too large (" + std::to_string(ldl.nnz_l()) + " nonzeros)";
        return out;
    }
    if (Clock::now() > deadline) {
        out.status = IpmStatus::TimeLimit;
        return out;
    }
    const double rho = kRho, delta = kDelta;
    Vec V(N + nqo + nza + M);
    for (int k = 0; k < nqo; ++k) V[N + k] = -P.qv[k];
    for (int p = 0; p < nza; ++p) V[N + nqo + p] = P.Ax[p];
    for (int i = 0; i < M; ++i) V[N + nqo + nza + i] = delta;
    Vec theta(N, 0);  // Theta^-1 = zl/xl + zu/xu
    auto factor = [&]() {
        for (int j = 0; j < N; ++j) V[j] = -(P.qd[j] + theta[j] + rho);
        if (std::getenv("IPM_DBG")) { double mx = 0, mn = INFINITY; for (int j = 0; j < N; ++j) mx = std::max(mx, theta[j]), mn = std::min(mn, theta[j]); std::fprintf(stderr, "  theta range %.3e %.3e\n", mn, mx); }
        return ldl.factor(V, deadline, true);
    };
    Vec t1, t2, t3;
    auto kkt_mul = [&](const Vec& v, Vec& r) {  // unregularised K v
        Vec vz(v.begin(), v.begin() + N), vy(v.begin() + N, v.end());
        P.mul_Q(vz, t1);
        P.mul_At(vy, t2);
        P.mul_A(vz, t3);
        r.resize(N + M);
        for (int j = 0; j < N; ++j) r[j] = -t1[j] - theta[j] * vz[j] + t2[j];
        for (int i = 0; i < M; ++i) r[N + i] = t3[i] + kDelta * kKeepDelta * vy[i];
    };
    Vec rr, dd, kv;
    double ref_res = 0;  // relative residual of the last refined solve (for the log)
    // Solve K sol = rhs with the unregularised K: restarted GMRES, right-preconditioned by the
    // regularised LDL^T factor (plain iterative refinement stalls when Theta^-1 << rho).
    const int gm = kRefMax;
    std::vector<Vec> Vb(gm + 1, Vec(N + M)), Zb(gm, Vec(N + M));
    std::vector<Vec> H(gm + 1, Vec(gm, 0));
    Vec cs(gm), sn(gm), gg(gm + 1), yy(gm);
    auto dot = [](const Vec& a, const Vec& b) {
        double s = 0;
        for (size_t k = 0; k < a.size(); ++k) s += a[k] * b[k];
        return s;
    };
    auto kkt_solve = [&](const Vec& rhs, Vec& sol) {
        const int n = N + M;
        sol = rhs;
        ldl.solve(sol);
        const double target = 1e-15 * (1 + norm_inf(rhs));
        double rn = INFINITY;
        for (int restart = 0; restart < 4; ++restart) {
            kkt_mul(sol, kv);
            rr.resize(n);
            for (int k = 0; k < n; ++k) rr[k] = rhs[k] - kv[k];
            rn = norm_inf(rr);
            if (!std::isfinite(rn) || rn <= target) break;
            const double beta = std::sqrt(dot(rr, rr));
            for (int k = 0; k < n; ++k) Vb[0][k] = rr[k] / beta;
            std::fill(gg.begin(), gg.end(), 0.0);
            gg[0] = beta;
            int used = 0;
            for (int j = 0; j < gm; ++j) {
                Zb[j] = Vb[j];
                ldl.solve(Zb[j]);
                kkt_mul(Zb[j], kv);
                for (int i = 0; i <= j; ++i) {  // modified Gram-Schmidt
                    double h = dot(kv, Vb[i]);
                    H[i][j] = h;
                    for (int k = 0; k < n; ++k) kv[k] -= h * Vb[i][k];
                }
                double hn = std::sqrt(dot(kv, kv));
                H[j + 1][j] = hn;
                for (int i = 0; i < j; ++i) {
                    double t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                    H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                    H[i][j] = t;
                }
                double den = std::hypot(H[j][j], H[j + 1][j]);
                if (!(den > 0)) break;
                cs[j] = H[j][j] / den, sn[j] = H[j + 1][j] / den;
                H[j][j] = den, H[j + 1][j] = 0;
                gg[j + 1] = -sn[j] * gg[j];
                gg[j] = cs[j] * gg[j];
                used = j + 1;
                if (!(hn > 0) || std::abs(gg[j + 1]) <= 0.1 * target) break;
                for (int k = 0; k < n; ++k) Vb[j + 1][k] = kv[k] / hn;
            }
            for (int i = used - 1; i >= 0; --i) {
                double t = gg[i];
                for (int k = i + 1; k < used; ++k) t -= H[i][k] * yy[k];
                yy[i] = t / H[i][i];
            }
            Vec trial = sol;
            for (int i = 0; i < used; ++i)
                for (int k = 0; k < n; ++k) trial[k] += yy[i] * Zb[i][k];
            kkt_mul(trial, kv);
            double tn = 0;
            for (int k = 0; k < n; ++k) tn = std::max(tn, std::abs(rhs[k] - kv[k]));
            if (!(tn < rn)) break;  // no improvement: keep the current solution
            sol.swap(trial);
            rn = tn;
            if (rn <= target || used < gm) break;
        }
        ref_res = rn / (1 + norm_inf(rhs));
    };

    const double bnorm = std::max({norm_inf(P.b), [&] {
                                       double m = 0;
                                       for (int j = 0; j < N; ++j) {
                                           if (hl[j]) m = std::max(m, std::abs(P.l[j]));
                                           if (hu[j]) m = std::max(m, std::abs(P.u[j]));
                                       }
                                       return m;
                                   }()});
    const double cnorm = std::max(norm_inf(P.c), norm_inf(P.qd));

    Vec z(N, 0), y(M, 0), xl(N, 0), xu(N, 0), zl(N, 0), zu(N, 0);
    Vec rhs(N + M), sol;
    // ---- starting point (Mehrotra, adapted to bounds) ----
    std::fill(theta.begin(), theta.end(), 1.0);
    if (!factor()) {
        out.status = IpmStatus::TimeLimit;
        return out;
    }
    for (int j = 0; j < N; ++j) rhs[j] = 0;
    for (int i = 0; i < M; ++i) rhs[N + i] = P.b[i];
    kkt_solve(rhs, sol);
    for (int j = 0; j < N; ++j) z[j] = sol[j];
    Vec Qz, Aty, Az;
    P.mul_Q(z, Qz);
    for (int j = 0; j < N; ++j) rhs[j] = P.c[j] + Qz[j];
    for (int i = 0; i < M; ++i) rhs[N + i] = 0;
    kkt_solve(rhs, sol);
    for (int i = 0; i < M; ++i) y[i] = sol[N + i];
    P.mul_At(y, Aty);
    {
        double minx = 0, minz = 0;
        for (int j = 0; j < N; ++j) {
            double s = P.c[j] + Qz[j] - Aty[j];
            if (hl[j]) xl[j] = z[j] - P.l[j], minx = std::min(minx, xl[j]);
            if (hu[j]) xu[j] = P.u[j] - z[j], minx = std::min(minx, xu[j]);
            if (hl[j] && hu[j]) zl[j] = std::max(s, 0.0), zu[j] = std::max(-s, 0.0);
            else if (hl[j]) zl[j] = s;
            else if (hu[j]) zu[j] = -s;
            if (hl[j]) minz = std::min(minz, zl[j]);
            if (hu[j]) minz = std::min(minz, zu[j]);
        }
        double dp = -1.5 * minx, dz = -1.5 * minz, xs = 0, sx = 0, sz = 0;
        for (int j = 0; j < N; ++j) {
            if (hl[j]) xl[j] += dp, zl[j] += dz, xs += xl[j] * zl[j], sx += xl[j], sz += zl[j];
            if (hu[j]) xu[j] += dp, zu[j] += dz, xs += xu[j] * zu[j], sx += xu[j], sz += zu[j];
        }
        double ex = 0, ez = 0;
        if (nc > 0) {
            if (sz <= 0 || sx <= 0 || xs <= 0) {
                ex = ez = 1;
            } else {
                ex = 0.5 * xs / sz;
                ez = 0.5 * xs / sx;
            }
        }
        for (int j = 0; j < N; ++j) {
            if (hl[j]) xl[j] = std::max(xl[j] + ex, 1e-8), zl[j] = std::max(zl[j] + ez, 1e-8);
            if (hu[j]) xu[j] = std::max(xu[j] + ex, 1e-8), zu[j] = std::max(zu[j] + ez, 1e-8);
        }
    }

    const bool lp = nqo == 0 && norm_inf(P.qd) == 0;
    Vec rp(M), rl(N), ru(N), rd(N), rcl(N), rcu(N);
    Vec dz(N), dy(M), dxl(N), dxu(N), dzl(N), dzu(N);
    Vec az(N), ay(M), axl(N), axu(N), azl(N), azu(N);  // affine direction
    auto newton = [&](Vec& Dz, Vec& Dy, Vec& Dxl, Vec& Dxu, Vec& Dzl, Vec& Dzu) {
        for (int j = 0; j < N; ++j) {
            double r1 = rd[j];
            if (hl[j]) r1 += (rcl[j] + zl[j] * rl[j]) / xl[j];
            if (hu[j]) r1 -= (rcu[j] - zu[j] * ru[j]) / xu[j];
            rhs[j] = -r1;
        }
        for (int i = 0; i < M; ++i) rhs[N + i] = rp[i];
        kkt_solve(rhs, sol);
        for (int j = 0; j < N; ++j) {
            Dz[j] = sol[j];
            Dxl[j] = hl[j] ? Dz[j] - rl[j] : 0;
            Dxu[j] = hu[j] ? ru[j] - Dz[j] : 0;
            Dzl[j] = hl[j] ? (rcl[j] - zl[j] * Dxl[j]) / xl[j] : 0;
            Dzu[j] = hu[j] ? (rcu[j] - zu[j] * Dxu[j]) / xu[j] : 0;
        }
        for (int i = 0; i < M; ++i) Dy[i] = sol[N + i];
    };
    auto steps = [&](const Vec& Dxl, const Vec& Dxu, const Vec& Dzl, const Vec& Dzu, double& ap, double& ad) {
        ap = 1, ad = 1;
        for (int j = 0; j < N; ++j) {
            if (hl[j]) {
                if (Dxl[j] < 0) ap = std::min(ap, -xl[j] / Dxl[j]);
                if (Dzl[j] < 0) ad = std::min(ad, -zl[j] / Dzl[j]);
            }
            if (hu[j]) {
                if (Dxu[j] < 0) ap = std::min(ap, -xu[j] / Dxu[j]);
                if (Dzu[j] < 0) ad = std::min(ad, -zu[j] / Dzu[j]);
            }
        }
    };

    // Certificate tests on a candidate (normalised internally).
    Vec g, Ad, Qd;
    auto farkas = [&](const Vec& yc) {  // y with b'y > sup_{l<=z<=u} (A'y)'z
        double ny = norm_inf(yc);
        if (!(ny > 0) || !std::isfinite(ny)) return false;
        P.mul_At(yc, g);
        double sup = 0, viol = 0, by = 0;
        for (int i = 0; i < M; ++i) by += P.b[i] * yc[i] / ny;
        for (int j = 0; j < N; ++j) {
            double gj = g[j] / ny;
            if (gj > 0) {
                if (hu[j]) sup += gj * P.u[j];
                else viol = std::max(viol, gj);
            } else if (gj < 0) {
                if (hl[j]) sup += gj * P.l[j];
                else viol = std::max(viol, -gj);
            }
        }
        return viol <= 1e-9 && by - sup > 1e-6 * (1 + 1e-3 * std::abs(by) + 1e-3 * std::abs(sup));
    };
    auto recession = [&](const Vec& dc) {  // A d = 0, Q d = 0, d in the bound cone, c'd < 0
        double nd = norm_inf(dc);
        if (!(nd > 0) || !std::isfinite(nd)) return false;
        Vec d(dc);
        for (double& v : d) v /= nd;
        double cd = 0, viol = 0;
        for (int j = 0; j < N; ++j) {
            cd += P.c[j] * d[j];
            if (hl[j] && d[j] < 0) viol = std::max(viol, -d[j]);
            if (hu[j] && d[j] > 0) viol = std::max(viol, d[j]);
        }
        if (viol > 1e-9 || cd > -1e-6 * (1 + cnorm)) return false;
        P.mul_A(d, Ad);
        P.mul_Q(d, Qd);
        if (!(norm_inf(Ad) <= 1e-9 && norm_inf(Qd) <= 1e-9)) return false;
        out.ray = d;
        return true;
    };

    int stall = 0, extra = 0;
    double last_ap = 0, last_ad = 0;
    double best_merit = INFINITY;
    for (long it = 0;; ++it) {
        out.iters = it;
        // residuals
        P.mul_A(z, Az);
        P.mul_Q(z, Qz);
        P.mul_At(y, Aty);
        double mu = 0, pobj = 0, dobj = 0;
        for (int i = 0; i < M; ++i) rp[i] = P.b[i] - Az[i], dobj += P.b[i] * y[i];
        for (int j = 0; j < N; ++j) {
            rl[j] = hl[j] ? P.l[j] - z[j] + xl[j] : 0;
            ru[j] = hu[j] ? P.u[j] - z[j] - xu[j] : 0;
            rd[j] = -(P.c[j] + Qz[j] - Aty[j] - zl[j] + zu[j]);
            if (hl[j]) mu += xl[j] * zl[j], dobj += P.l[j] * zl[j];
            if (hu[j]) mu += xu[j] * zu[j], dobj -= P.u[j] * zu[j];
            pobj += P.c[j] * z[j] + 0.5 * z[j] * Qz[j];
            dobj -= 0.5 * z[j] * Qz[j];
        }
        mu = nc ? mu / nc : 0;
        const double pres = std::max({norm_inf(rp), norm_inf(rl), norm_inf(ru)}) / (1 + bnorm);
        const double dres = norm_inf(rd) / (1 + cnorm);
        if (std::getenv("IPM_DBG")) std::fprintf(stderr, "  rp %.2e rl %.2e ru %.2e rd %.2e y %.2e z %.2e zl %.2e zu %.2e\n", norm_inf(rp), norm_inf(rl), norm_inf(ru), norm_inf(rd), norm_inf(y), norm_inf(z), norm_inf(zl), norm_inf(zu));
        if (std::getenv("IPM_DBG") && it % 10 == 9) { int jm = 0; for (int j = 0; j < N; ++j) if (std::abs(rd[j]) > std::abs(rd[jm])) jm = j; int nzc = P.Ap[jm + 1] - P.Ap[jm]; std::fprintf(stderr, "  maxrd j %d rd %.3e z %.3e l %.3e u %.3e xl %.3e zl %.3e xu %.3e zu %.3e nnzcol %d c %.3e\n", jm, rd[jm], z[jm], P.l[jm], P.u[jm], xl[jm], zl[jm], xu[jm], zu[jm], nzc, P.c[jm]); }
        const double gap = std::abs(pobj - dobj) / (1 + std::abs(pobj));
        if (opt.verbose)
            std::fprintf(stderr, "it %3ld pobj %+.10e dobj %+.10e pres %.2e dres %.2e gap %.2e mu %.2e reg %d"
                         " ap %.2e ad %.2e kkt %.1e\n", it, pobj, dobj, pres, dres, gap, mu, ldl.nreg, last_ap,
                         last_ad, ref_res);
        if (!std::isfinite(pres + dres + gap + mu)) {
            out.msg = "non-finite iterate";
            break;
        }
        out.z = z, out.y = y;
        if (std::max({pres, dres, gap}) <= 10 * opt.tol) {
            if (accept(z, y)) {
                out.status = IpmStatus::Optimal;
                return out;
            }
            if (std::max({pres, dres, gap}) <= opt.tol && ++extra > 15) {
                out.msg = "converged on the scaled problem but not on the original model";
                break;
            }
        }
        // certificates
        if (farkas(y) || (it > 0 && farkas(dy))) {
            out.status = IpmStatus::Infeasible;
            out.msg = "Farkas certificate";
            return out;
        }
        if (recession(z) || (it > 0 && recession(dz))) {
            out.status = IpmStatus::DualInfeasible;  // ipm_solve promotes it to Unbounded only after re-checking on the original model
            out.msg = "recession direction";
            return out;
        }
        if (it >= opt.max_iter) {
            out.status = IpmStatus::IterationLimit;
            break;
        }
        if (Clock::now() > deadline) {
            out.status = IpmStatus::TimeLimit;
            break;
        }
        // Newton step
        for (int j = 0; j < N; ++j) theta[j] = (hl[j] ? zl[j] / xl[j] : 0) + (hu[j] ? zu[j] / xu[j] : 0);
        if (!factor()) {
            out.status = IpmStatus::TimeLimit;
            break;
        }
        for (int j = 0; j < N; ++j) rcl[j] = hl[j] ? -xl[j] * zl[j] : 0, rcu[j] = hu[j] ? -xu[j] * zu[j] : 0;
        newton(az, ay, axl, axu, azl, azu);
        double ap, ad;
        steps(axl, axu, azl, azu, ap, ad);
        if (lp == false) ap = ad = std::min(ap, ad);
        if (nc > 0) {
            double mua = 0;
            for (int j = 0; j < N; ++j) {
                if (hl[j]) mua += (xl[j] + ap * axl[j]) * (zl[j] + ad * azl[j]);
                if (hu[j]) mua += (xu[j] + ap * axu[j]) * (zu[j] + ad * azu[j]);
            }
            mua /= nc;
            double sigma = std::pow(std::max(0.0, mua) / mu, 3);
            sigma = std::min(sigma, 1.0);
            for (int j = 0; j < N; ++j) {
                if (hl[j]) rcl[j] = sigma * mu - xl[j] * zl[j] - axl[j] * azl[j];
                if (hu[j]) rcu[j] = sigma * mu - xu[j] * zu[j] - axu[j] * azu[j];
            }
            newton(dz, dy, dxl, dxu, dzl, dzu);
            steps(dxl, dxu, dzl, dzu, ap, ad);
        } else {
            dz = az, dy = ay;
            ap = ad = 1;
        }
        ap = std::min(1.0, 0.995 * ap);
        ad = std::min(1.0, 0.995 * ad);
        if (nc == 0) ap = ad = 1;
        if (!lp) ap = ad = std::min(ap, ad);
        for (int j = 0; j < N; ++j) {
            z[j] += ap * dz[j];
            if (hl[j]) xl[j] += ap * dxl[j], zl[j] += ad * dzl[j];
            if (hu[j]) xu[j] += ap * dxu[j], zu[j] += ad * dzu[j];
        }
        for (int i = 0; i < M; ++i) y[i] += ad * dy[i];
        last_ap = ap, last_ad = ad;
        double merit = std::max({pres, dres, gap});
        if (merit < 0.9 * best_merit) best_merit = merit, stall = 0;
        else if (++stall > 30) {
            out.msg = "no progress in 30 iterations";
            break;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------------------------
struct Measures {
    double pres = 0, dres = 0, gap = 0, pobj = 0, dobj = 0, row_viol = 0, bound_viol = 0;
    Vec act, rc;
};

// Original-model KKT measures in minimisation form (c, Q already multiplied by the sense).
Measures measure(const Model& md, const Vec& c, const std::vector<QEntry>& Q, const Vec& clo, const Vec& cup,
                 const Vec& rlo, const Vec& rup, const Vec& x, const Vec& y) {
    const size_t n = x.size(), m = rlo.size();
    Measures r;
    r.act.assign(m, 0);
    for (size_t j = 0; j < n; ++j)
        for (const Entry& e : md.cols[j]) r.act[e.index] += e.value * x[j];
    auto viol = [](double v, double lo, double up, double& absv) {
        double a = std::max({lo - v, v - up, 0.0});
        absv = std::max(absv, a);
        if (a == 0) return 0.0;
        return a / (1 + std::abs(v < lo ? lo : up));
    };
    for (size_t i = 0; i < m; ++i) r.pres = std::max(r.pres, viol(r.act[i], rlo[i], rup[i], r.row_viol));
    for (size_t j = 0; j < n; ++j) r.pres = std::max(r.pres, viol(x[j], clo[j], cup[j], r.bound_viol));
    Vec Qx(n, 0);
    for (const QEntry& q : Q) {
        Qx[q.row] += q.value * x[q.col];
        if (q.row != q.col) Qx[q.col] += q.value * x[q.row];
    }
    r.rc.assign(n, 0);
    double xQx = 0, dv = 0;
    for (size_t j = 0; j < n; ++j) {
        double s = c[j] + Qx[j];
        for (const Entry& e : md.cols[j]) s -= e.value * y[e.index];
        r.rc[j] = s;
        xQx += x[j] * Qx[j];
        r.pobj += c[j] * x[j];
    }
    r.pobj += 0.5 * xQx;
    r.dobj = -0.5 * xQx;
    // A multiplier v on a "lo <= t <= up" constraint must be >= 0 if only lo is finite, <= 0 if only
    // up is finite, 0 if neither; its dual objective term is v*lo (v > 0) or v*up (v < 0).
    auto dual_term = [&](double v, double lo, double up) {
        if (v > 0) {
            if (std::isfinite(lo)) r.dobj += v * lo;
            else dv = std::max(dv, v);
        } else if (v < 0) {
            if (std::isfinite(up)) r.dobj += v * up;
            else dv = std::max(dv, -v);
        }
    };
    for (size_t i = 0; i < m; ++i) dual_term(y[i], rlo[i], rup[i]);
    for (size_t j = 0; j < n; ++j) dual_term(r.rc[j], clo[j], cup[j]);
    r.dres = dv / (1 + norm_inf(c));
    r.pobj += md.obj_const * (md.maximize ? -1 : 1);
    r.dobj += md.obj_const * (md.maximize ? -1 : 1);
    r.gap = std::abs(r.pobj - r.dobj) / (1 + std::abs(r.pobj));
    return r;
}

// Inertia test: Q + tau I with tau = 1e-7 max|Q| must have only positive pivots.
bool q_convex(int n, const std::vector<QEntry>& q) {
    if (q.empty()) return true;
    std::vector<int> id(n, -1);
    int k = 0;
    double qmax = 0;
    for (const QEntry& e : q) {
        for (int v : {e.row, e.col})
            if (id[v] < 0) id[v] = k++;
        qmax = std::max(qmax, std::abs(e.value));
    }
    std::vector<int> ei, ej, es;
    Vec V(k + q.size(), 0);
    for (int i = 0; i < k; ++i) ei.push_back(i), ej.push_back(i), es.push_back(i), V[i] = 1e-7 * qmax;
    for (size_t t = 0; t < q.size(); ++t) {
        if (q[t].row == q[t].col) {
            V[id[q[t].row]] += q[t].value;
        } else {
            ei.push_back(id[q[t].row]), ej.push_back(id[q[t].col]), es.push_back(k + static_cast<int>(t));
            V[k + t] = q[t].value;
        }
    }
    Ldl f;
    f.analyze(k, ei, ej, es, std::vector<signed char>(k, 1));
    f.factor(V, Clock::time_point::max(), false);
    return f.nwrong == 0;
}

}  // namespace

const char* ipm_status_name(IpmStatus s) {
    switch (s) {
        case IpmStatus::Optimal: return "optimal";
        case IpmStatus::Infeasible: return "infeasible";
        case IpmStatus::Unbounded: return "unbounded";
        case IpmStatus::DualInfeasible: return "dual_infeasible";
        case IpmStatus::TimeLimit: return "time_limit";
        case IpmStatus::IterationLimit: return "iteration_limit";
        case IpmStatus::NumericalFailure: return "numerical_failure";
        case IpmStatus::Nonconvex: return "nonconvex";
        case IpmStatus::Unsupported: return "unsupported";
    }
    return "?";
}

IpmResult ipm_solve(const Model& md, const IpmOptions& opt) {
    const auto t0 = Clock::now();
    const auto deadline = t0 + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(opt.time_limit));
    IpmResult res;
    const int n0 = static_cast<int>(md.col_names.size()), m0 = static_cast<int>(md.row_names.size());
    if (md.has_integers()) {
        res.status = IpmStatus::Unsupported;
        res.message = "integer columns are not supported by the IPM";
        return res;
    }
    const double sense = md.maximize ? -1 : 1;
    auto lo_of = [](double v) { return v <= -kBig ? -INFINITY : v; };
    auto up_of = [](double v) { return v >= kBig ? INFINITY : v; };
    Vec clo(n0), cup(n0), rlo(m0), rup(m0), c(n0);
    for (int j = 0; j < n0; ++j) clo[j] = lo_of(md.col_lo[j]), cup[j] = up_of(md.col_up[j]), c[j] = sense * md.cost[j];
    for (int i = 0; i < m0; ++i) rlo[i] = lo_of(md.row_lo[i]), rup[i] = up_of(md.row_up[i]);
    std::vector<QEntry> Q = md.qobj;
    for (QEntry& q : Q) q.value *= sense;
    res.x.assign(n0, 0);
    res.row_dual.assign(m0, 0);

    auto finish = [&](const Vec& x, const Vec& y_min) {  // fill measures and model-sense outputs
        Measures ms = measure(md, c, Q, clo, cup, rlo, rup, x, y_min);
        res.x = x;
        res.row_activity = ms.act;
        res.row_dual.resize(m0);
        res.reduced_cost.resize(n0);
        for (int i = 0; i < m0; ++i) res.row_dual[i] = sense * y_min[i];
        for (int j = 0; j < n0; ++j) res.reduced_cost[j] = sense * ms.rc[j];
        res.objective = sense * ms.pobj;
        res.dual_objective = sense * ms.dobj;
        res.primal_res = ms.pres, res.dual_res = ms.dres, res.gap = ms.gap;
        res.max_row_viol = ms.row_viol, res.max_bound_viol = ms.bound_viol;
        return ms;
    };
    auto infeasible = [&](const std::string& why) {
        res.status = IpmStatus::Infeasible;
        res.message = why;
        finish(res.x, res.row_dual);
        return res;
    };
    for (int j = 0; j < n0; ++j)
        if (clo[j] > cup[j] || clo[j] == INFINITY || cup[j] == -INFINITY)
            return infeasible("column " + md.col_names[j] + " has empty bounds");
    for (int i = 0; i < m0; ++i)
        if (rlo[i] > rup[i] || rlo[i] == INFINITY || rup[i] == -INFINITY)
            return infeasible("row " + md.row_names[i] + " has empty bounds");

    // ---- light presolve ----
    std::vector<char> fixed(n0);
    Vec x0(n0, 0), cadj = c;
    for (int j = 0; j < n0; ++j)
        if (clo[j] == cup[j]) fixed[j] = 1, x0[j] = clo[j];
    for (const QEntry& q : Q) {
        if (fixed[q.row] && !fixed[q.col]) cadj[q.col] += q.value * x0[q.row];
        else if (fixed[q.col] && !fixed[q.row]) cadj[q.row] += q.value * x0[q.col];
    }
    Vec off(m0, 0);
    std::vector<int> cnt(m0, 0);
    for (int j = 0; j < n0; ++j)
        for (const Entry& e : md.cols[j]) {
            if (fixed[j]) off[e.index] += e.value * x0[j];
            else ++cnt[e.index];
        }
    std::vector<int> rmap(m0, -1), rows;
    for (int i = 0; i < m0; ++i) {
        if (cnt[i] == 0) {
            double tol = 1e-9 * (1 + std::abs(off[i]));
            if (off[i] < rlo[i] - tol || off[i] > rup[i] + tol)
                return infeasible("row " + md.row_names[i] + " has no free columns and is violated");
            continue;
        }
        if (!std::isfinite(rlo[i]) && !std::isfinite(rup[i])) continue;
        rmap[i] = static_cast<int>(rows.size());
        rows.push_back(i);
    }
    std::vector<char> qtouch(n0, 0);
    for (const QEntry& q : Q)
        if (!fixed[q.row] && !fixed[q.col]) qtouch[q.row] = qtouch[q.col] = 1;
    std::vector<int> cmap(n0, -1), cols;
    std::string ray_col;  // an empty column whose cost pushes it to an infinite bound
    for (int j = 0; j < n0; ++j) {
        if (fixed[j]) continue;
        bool empty = !qtouch[j];
        for (const Entry& e : md.cols[j]) empty = empty && rmap[e.index] < 0;
        if (!empty) {
            cmap[j] = static_cast<int>(cols.size());
            cols.push_back(j);
            continue;
        }
        double v;
        if (cadj[j] > 0) v = clo[j];
        else if (cadj[j] < 0) v = cup[j];
        else v = std::min(std::max(0.0, clo[j]), cup[j]);
        if (!std::isfinite(v)) {
            ray_col = md.col_names[j];
            v = std::isfinite(clo[j]) ? clo[j] : std::isfinite(cup[j]) ? cup[j] : 0;
        }
        x0[j] = v;
    }
    const int nk = static_cast<int>(cols.size()), mk = static_cast<int>(rows.size());
    std::vector<QEntry> Qk;
    for (const QEntry& q : Q)
        if (cmap[q.row] >= 0 && cmap[q.col] >= 0) Qk.push_back({cmap[q.row], cmap[q.col], q.value});
    if (!q_convex(nk, Qk)) {
        res.status = IpmStatus::Nonconvex;
        res.message = "Q is not positive semidefinite in the minimisation sense";
        return res;
    }
    // structural CSC (kept rows x kept cols)
    std::vector<int> Kp(nk + 1, 0), Ki;
    Vec Kx;
    for (int k = 0; k < nk; ++k) {
        for (const Entry& e : md.cols[cols[k]])
            if (rmap[e.index] >= 0) Ki.push_back(rmap[e.index]), Kx.push_back(e.value);
        Kp[k + 1] = static_cast<int>(Ki.size());
    }
    // ---- scaling: geometric passes on A, then Ruiz equilibration of [Q A'; A 0], powers of two ----
    Vec R(mk, 1), D(nk, 1);
    for (int pass = 0; pass < 6; ++pass) {
        for (int k = 0; k < nk; ++k) {
            double mn = INFINITY, mx = 0;
            for (int p = Kp[k]; p < Kp[k + 1]; ++p) {
                double a = std::abs(Kx[p]) * R[Ki[p]];
                mn = std::min(mn, a), mx = std::max(mx, a);
            }
            if (mx > 0) D[k] = 1 / std::sqrt(mn * mx);
        }
        Vec mn(mk, INFINITY), mx(mk, 0);
        for (int k = 0; k < nk; ++k)
            for (int p = Kp[k]; p < Kp[k + 1]; ++p) {
                double a = std::abs(Kx[p]) * D[k];
                mn[Ki[p]] = std::min(mn[Ki[p]], a), mx[Ki[p]] = std::max(mx[Ki[p]], a);
            }
        for (int i = 0; i < mk; ++i)
            if (mx[i] > 0) R[i] = 1 / std::sqrt(mn[i] * mx[i]);
    }
    for (int pass = 0; pass < 20; ++pass) {
        Vec cn(nk, 0), rn(mk, 0);
        for (int k = 0; k < nk; ++k)
            for (int p = Kp[k]; p < Kp[k + 1]; ++p) {
                double a = std::abs(Kx[p]) * R[Ki[p]] * D[k];
                cn[k] = std::max(cn[k], a), rn[Ki[p]] = std::max(rn[Ki[p]], a);
            }
        for (const QEntry& q : Qk) {
            double a = std::abs(q.value) * D[q.row] * D[q.col];
            cn[q.row] = std::max(cn[q.row], a), cn[q.col] = std::max(cn[q.col], a);
        }
        for (int k = 0; k < nk; ++k)
            if (cn[k] > 0) D[k] /= std::sqrt(cn[k]);
        for (int i = 0; i < mk; ++i)
            if (rn[i] > 0) R[i] /= std::sqrt(rn[i]);
    }
    for (double& d : D) d = std::exp2(std::round(std::log2(d)));
    for (double& r : R) r = std::exp2(std::round(std::log2(r)));
    double cmax = 0;
    for (int k = 0; k < nk; ++k) cmax = std::max(cmax, std::abs(cadj[cols[k]]) * D[k]);
    for (const QEntry& q : Qk) cmax = std::max(cmax, std::abs(q.value) * D[q.row] * D[q.col]);
    const double sig = std::exp2(-std::round(std::log2(std::max(1.0, cmax))));

    // ---- internal problem ----
    Qp P;
    std::vector<int> slack_row;  // internal slack -> kept row
    for (int i = 0; i < mk; ++i)
        if (rlo[rows[i]] != rup[rows[i]]) slack_row.push_back(i);
    P.N = nk + static_cast<int>(slack_row.size());
    P.M = mk;
    P.Ap = Kp;
    P.Ai = Ki;
    P.Ax.resize(Kx.size());
    for (int k = 0; k < nk; ++k)
        for (int p = Kp[k]; p < Kp[k + 1]; ++p) P.Ax[p] = Kx[p] * R[Ki[p]] * D[k];
    for (int i : slack_row) P.Ai.push_back(i), P.Ax.push_back(-1), P.Ap.push_back(static_cast<int>(P.Ai.size()));
    P.qd.assign(P.N, 0);
    for (const QEntry& q : Qk) {
        double v = sig * q.value * D[q.row] * D[q.col];
        if (q.row == q.col) P.qd[q.row] += v;
        else P.qr.push_back(q.row), P.qc.push_back(q.col), P.qv.push_back(v);
    }
    P.c.assign(P.N, 0);
    P.l.resize(P.N), P.u.resize(P.N);
    for (int k = 0; k < nk; ++k) {
        int j = cols[k];
        P.c[k] = sig * cadj[j] * D[k];
        P.l[k] = clo[j] / D[k], P.u[k] = cup[j] / D[k];
    }
    P.b.assign(mk, 0);
    for (int i = 0; i < mk; ++i)
        if (rlo[rows[i]] == rup[rows[i]]) P.b[i] = R[i] * (rlo[rows[i]] - off[rows[i]]);
    for (size_t s = 0; s < slack_row.size(); ++s) {
        int i = slack_row[s], r = rows[i];
        P.l[nk + s] = R[i] * (rlo[r] - off[r]);
        P.u[nk + s] = R[i] * (rup[r] - off[r]);
    }
    res.kkt_dim = P.N + P.M;

    // map an internal (z, y) to the original model
    auto to_orig = [&](const Vec& z, const Vec& y, Vec& x, Vec& yo) {
        x = x0;
        for (int k = 0; k < nk; ++k) x[cols[k]] = std::min(std::max(z[k] * D[k], clo[cols[k]]), cup[cols[k]]);
        yo.assign(m0, 0);
        for (int i = 0; i < mk; ++i) yo[rows[i]] = R[i] * y[i] / sig;
    };
    Vec x, yo;
    auto accept = [&](const Vec& z, const Vec& y) {
        to_orig(z, y, x, yo);
        Measures ms = measure(md, c, Q, clo, cup, rlo, rup, x, yo);
        // With a ray column the original model has no dual solution (that column is dual infeasible
        // by construction), so only primal feasibility can be confirmed here; the caller reports unbounded.
        if (!ray_col.empty()) return ms.pres <= opt.tol;
        return ms.pres <= opt.tol && ms.dres <= opt.tol && ms.gap <= opt.tol;
    };
    RunOut ro;
    if (P.N > 0) {
        ro = run_ipm(P, opt, deadline, accept);
    } else {
        ro.status = IpmStatus::Optimal;
    }
    res.iterations = ro.iters;
    res.factor_nnz = ro.fnnz;
    res.message = ro.msg;
    if (ro.z.size() == static_cast<size_t>(P.N) && ro.y.size() == static_cast<size_t>(P.M)) to_orig(ro.z, ro.y, x, yo);
    else x = x0, yo.assign(m0, 0);
    Measures ms = finish(x, yo);
    res.status = ro.status;
    if (res.status == IpmStatus::DualInfeasible) {
        // "Unbounded" needs two facts, both checked here on the original unscaled model and nowhere else:
        //   (1) a feasible point (row and column violation <= opt.tol, the measure used for Optimal), and
        //   (2) an improving ray d (max-norm 1): d in the column cone (0 for a boxed column, >= 0 with only a lower bound,
        //       <= 0 with only an upper bound), A d in the row cone (0 for a two-sided row, >= 0 / <= 0 for one-sided),
        //       Q d = 0 and c'd < -1e-6 (1 + |c|inf), each with slack opt.tol (row slack scaled by the row's sum |a_j d_j|).
        // If the returned iterate is not feasible, a zero-objective solve of the same constraints (c = 0, Q = 0) supplies a
        // point; it must pass the same check. Anything that fails stays DualInfeasible.
        auto ray_ok = [&]() {
            const double rtol = opt.tol;  // same tolerance as the point check
            if (ro.ray.size() != static_cast<size_t>(P.N) || nk == 0) return false;  // first nk entries: columns, rest: slacks
            Vec d(n0, 0);
            double nd = 0;
            for (int k = 0; k < nk; ++k) d[cols[k]] = ro.ray[k] * D[k], nd = std::max(nd, std::abs(d[cols[k]]));
            if (!(nd > 0) || !std::isfinite(nd)) return false;
            for (double& v : d) v /= nd;
            for (int j = 0; j < n0; ++j) {
                bool lo = std::isfinite(clo[j]), up = std::isfinite(cup[j]);
                if ((lo && d[j] < -rtol) || (up && d[j] > rtol)) return false;
                if (lo && up && std::abs(d[j]) > rtol) return false;
            }
            Vec ad(m0, 0), sc(m0, 0);
            for (int j = 0; j < n0; ++j)
                for (const Entry& e : md.cols[j]) ad[e.index] += e.value * d[j], sc[e.index] += std::abs(e.value * d[j]);
            for (int i = 0; i < m0; ++i) {
                bool lo = std::isfinite(rlo[i]), up = std::isfinite(rup[i]);
                double tol = rtol * (1 + sc[i]);
                if ((lo && ad[i] < -tol) || (up && ad[i] > tol)) return false;
                if (lo && up && std::abs(ad[i]) > tol) return false;
            }
            Vec qd(n0, 0), qs(n0, 0);
            for (const QEntry& q : Q) {
                qd[q.row] += q.value * d[q.col], qs[q.row] += std::abs(q.value * d[q.col]);
                if (q.row != q.col) qd[q.col] += q.value * d[q.row], qs[q.col] += std::abs(q.value * d[q.row]);
            }
            for (int j = 0; j < n0; ++j)
                if (std::abs(qd[j]) > rtol * (1 + qs[j])) return false;
            double cd = 0;
            for (int j = 0; j < n0; ++j) cd += c[j] * d[j];
            return cd < -1e-6 * (1 + norm_inf(c));
        };
        if (ray_ok()) {
            if (ms.pres <= opt.tol) {
                res.status = IpmStatus::Unbounded;
                res.message = "recession direction and returned point verified on the original model";
            } else {
                Model fm = md;  // same constraints, no objective
                std::fill(fm.cost.begin(), fm.cost.end(), 0.0);
                fm.qobj.clear();
                fm.obj_const = 0;
                fm.maximize = false;
                IpmOptions fo = opt;
                fo.time_limit = std::max(0.0, std::chrono::duration<double>(deadline - Clock::now()).count());
                IpmResult fr = ipm_solve(fm, fo);
                if (fr.status == IpmStatus::Optimal && fr.x.size() == static_cast<size_t>(n0) &&
                    measure(md, c, Q, clo, cup, rlo, rup, fr.x, Vec(m0, 0)).pres <= opt.tol) {
                    res.status = IpmStatus::Unbounded;
                    res.message = "recession direction verified; feasible point from a zero-objective solve, verified on the original model";
                    finish(fr.x, Vec(m0, 0));
                }
            }
        }
    }
    if (res.status == IpmStatus::Optimal && !ray_col.empty()) {
        // The ray column is dual infeasible by construction, so only primal feasibility of the rest is
        // checkable: feasible point plus a free improving ray means unbounded.
        if (ms.pres <= opt.tol) {
            res.status = IpmStatus::Unbounded;
            res.message = "empty column " + ray_col + " is unbounded in its cost direction";
        } else {
            res.status = IpmStatus::NumericalFailure;
            res.message = "feasibility of the reduced problem could not be confirmed";
        }
    } else if (res.status == IpmStatus::Optimal) {
        if (!(ms.pres <= opt.tol && ms.dres <= opt.tol && ms.gap <= opt.tol)) {  // e.g. the P.N == 0 path
            res.status = IpmStatus::NumericalFailure;
            res.message = "final original-model recheck failed";
        }
    }
    (void)t0;
    return res;
}
