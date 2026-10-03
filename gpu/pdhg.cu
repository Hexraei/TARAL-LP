// TARAL-LP GPU track: PDHG (primal-dual hybrid gradient, PDLP-style restarted average) LP solver. PROTOTYPE.
// Near-optimal to a tolerance, NOT an exact solve. CUDA runtime + C++17 standard library only: own kernels for CSR
// matvec (A and A^T resident on the device), projections, averages and reductions; no cuSPARSE/cuBLAS/Thrust/CUB.
// The same algorithm runs on the CPU (--device cpu, std::thread pool) so iterations and times compare like for like.
//
// LP: min c'x + obj_const  s.t.  row_lo <= A x <= row_up,  col_lo <= x <= col_up   (src/taral.hpp Model semantics).
// Algorithm: Ruiz (10 passes, inf-norm) + Pock-Chambolle (alpha 1) scaling; constant step eta = 0.95/||A||_2 from
// power iteration; tau = eta/w, sigma = eta*w with primal weight w updated at restarts (smoothing 0.5); uniform
// averages; KKT-based adaptive restarts to the better of current/average (beta 0.2 / 0.8 / artificial 0.36).
// Every --eval-freq iterations the KKT error of current and average is reduced on the device to a few scalars.
// Stop when, on the ORIGINAL (unscaled) problem, ||primal res||_2 <= tol (1 + ||b||_2), ||dual res||_2 <= tol (1 +
// ||c||_2) and |pobj - dobj| <= tol (1 + |pobj| + |dobj|), and, only when --row-rel K is given (default 0 = off), the per-row relative primal
// residual ||viol_i / (1 + |b_i|)||_2 <= K tol, so one row cannot hide behind a large ||b|| (pds-100, fome21).
//
// Build: nvcc -O3 -arch=sm_120 -std=c++17 -Xcompiler "-O3 -march=native -pthread" -o out/pdhg gpu/pdhg.cu src/mps.cpp
// Run:   out/pdhg MODEL.mps --device gpu|cpu [--threads N] [--tol 1e-4] [--time-limit S] [--max-iter N]
//                 [--json OUT.json] [--sol OUT.sol] [--dual OUT.dual] [--fp32] [--eval-freq 64] [--row-rel 0] [--no-graph] [--verbose]
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <condition_variable>
#include <mutex>
#include <functional>
#include <string>
#include <thread>
#include <vector>

#include "../src/taral.hpp"

#define CHECK(call)                                                                                     \
    do {                                                                                                \
        cudaError_t e_ = (call);                                                                        \
        if (e_ != cudaSuccess) {                                                                        \
            std::fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
            std::exit(1);                                                                               \
        }                                                                                               \
    } while (0)

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// ---------------------------------------------------------------- host problem + scaling (shared by both backends)
struct Csr {
    int rows = 0;
    std::vector<int> ptr, idx;
    std::vector<double> val;
};

static Csr transpose(const Csr& a, int cols) {
    Csr t;
    t.rows = cols;
    t.ptr.assign(cols + 1, 0);
    for (int k : a.idx) t.ptr[k + 1]++;
    for (int j = 0; j < cols; ++j) t.ptr[j + 1] += t.ptr[j];
    t.idx.resize(a.idx.size());
    t.val.resize(a.val.size());
    std::vector<int> pos(t.ptr.begin(), t.ptr.end() - 1);
    for (int i = 0; i < a.rows; ++i)
        for (int k = a.ptr[i]; k < a.ptr[i + 1]; ++k) {
            int p = pos[a.idx[k]]++;
            t.idx[p] = i, t.val[p] = a.val[k];
        }
    return t;
}

// Scaled problem: A_s = diag(Dr) A diag(Dc), x = Dc x_s, A x = A_s x_s / Dr, z = Dr z_s, reduced cost d = d_s / Dc.
struct Scaled {
    int m = 0, n = 0;
    Csr A, AT;
    std::vector<double> c, l, u, rl, ru, Dr, Dc;
    double bnorm = 0, cnorm = 0;  // ORIGINAL-space norms used by the relative stopping test
};

static Scaled build_scaled(const Model& md, double sign) {
    Scaled P;
    P.n = int(md.cols.size()), P.m = int(md.row_lo.size());
    Csr T;  // A^T straight from the column-wise model
    T.rows = P.n;
    T.ptr.assign(P.n + 1, 0);
    for (int j = 0; j < P.n; ++j) {
        for (const Entry& e : md.cols[j])
            if (e.value != 0) T.idx.push_back(e.index), T.val.push_back(e.value);
        T.ptr[j + 1] = int(T.idx.size());
    }
    P.A = transpose(T, P.m);
    Csr& A = P.A;
    P.Dr.assign(P.m, 1), P.Dc.assign(P.n, 1);
    std::vector<double> rs(P.m), cs(P.n);
    auto apply = [&] {  // divide by sqrt of the row and column measures, accumulate the factors
        for (double& v : rs) v = v > 0 ? std::sqrt(v) : 1;
        for (double& v : cs) v = v > 0 ? std::sqrt(v) : 1;
        for (int i = 0; i < P.m; ++i)
            for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) A.val[k] /= rs[i] * cs[A.idx[k]];
        for (int i = 0; i < P.m; ++i) P.Dr[i] /= rs[i];
        for (int j = 0; j < P.n; ++j) P.Dc[j] /= cs[j];
    };
    for (int pass = 0; pass < 10; ++pass) {  // Ruiz equilibration, inf-norm
        std::fill(rs.begin(), rs.end(), 0), std::fill(cs.begin(), cs.end(), 0);
        for (int i = 0; i < P.m; ++i)
            for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) {
                double a = std::fabs(A.val[k]);
                rs[i] = std::max(rs[i], a), cs[A.idx[k]] = std::max(cs[A.idx[k]], a);
            }
        apply();
    }
    std::fill(rs.begin(), rs.end(), 0), std::fill(cs.begin(), cs.end(), 0);  // Pock-Chambolle, alpha = 1
    for (int i = 0; i < P.m; ++i)
        for (int k = A.ptr[i]; k < A.ptr[i + 1]; ++k) rs[i] += std::fabs(A.val[k]), cs[A.idx[k]] += std::fabs(A.val[k]);
    apply();
    P.AT = transpose(A, P.n);
    P.c.resize(P.n), P.l.resize(P.n), P.u.resize(P.n), P.rl.resize(P.m), P.ru.resize(P.m);
    for (int j = 0; j < P.n; ++j) {
        double c = sign * md.cost[j];
        P.cnorm += c * c;
        P.c[j] = c * P.Dc[j], P.l[j] = md.col_lo[j] / P.Dc[j], P.u[j] = md.col_up[j] / P.Dc[j];
    }
    for (int i = 0; i < P.m; ++i) {
        double lo = md.row_lo[i], up = md.row_up[i], b = 0;
        if (std::isfinite(lo) && std::isfinite(up)) b = std::max(std::fabs(lo), std::fabs(up));
        else if (std::isfinite(lo)) b = std::fabs(lo);
        else if (std::isfinite(up)) b = std::fabs(up);
        P.bnorm += b * b;
        P.rl[i] = lo * P.Dr[i], P.ru[i] = up * P.Dr[i];
    }
    P.bnorm = std::sqrt(P.bnorm), P.cnorm = std::sqrt(P.cnorm);
    return P;
}

// Scalars of one KKT evaluation (row part + column part), s = scaled space, o = original space.
enum { PR2S, PR2O, PR2R, DROW, DZ2, NROW };
// 1/(1+|b_i|) in original row units (b_i = larger finite |row bound|): weight of row i in the per-row-relative residual PR2R.
template <class T>
__host__ __device__ inline double row_rel_weight(double lo, double up, T dr) {
    double b = 0;
    if (lo > -INFINITY) b = lo < 0 ? -lo : lo;
    if (up < INFINITY) { double a = up < 0 ? -up : up; b = a > b ? a : b; }
    return 1.0 / (1.0 + b / double(dr));
}
enum { DR2S, DR2O, POBJ, DCOL, DX2, NCOL };
struct Kkt {
    double r[NROW] = {}, c[NCOL] = {};
    double pobj() const { return c[POBJ]; }
    double dobj() const { return r[DROW] + c[DCOL]; }
};

// Shared deterministic power-iteration start vector.
static double start_vec(int j) { return 0.5 + double((uint32_t(j) * 2654435761u) % 1000u) / 1000.0; }

// ---------------------------------------------------------------- GPU backend
template <class T>
__device__ inline T clampv(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }

template <int G, class T>
__device__ inline T group_dot(bool ok, int row, int lane, const int* __restrict__ ptr, const int* __restrict__ idx,
                              const T* __restrict__ val, const T* __restrict__ x) {
    T s = 0;
    if (ok)
        for (int k = ptr[row] + lane; k < ptr[row + 1]; k += G) s += val[k] * x[idx[k]];
    for (int off = G / 2; off > 0; off /= 2) s += __shfl_down_sync(0xffffffffu, s, off, G);
    return s;
}

// y = A x with G threads per row (G = 32 is the warp-per-row kernel of gpu/csr_matvec.cu).
template <int G, class T>
__global__ void k_spmv(int rows, const int* ptr, const int* idx, const T* val, const T* x, T* y) {
    int t = blockIdx.x * blockDim.x + threadIdx.x, row = t / G, lane = t % G;
    bool ok = row < rows;
    T s = group_dot<G>(ok, row, lane, ptr, idx, val, x);
    if (ok && lane == 0) y[row] = s;
}

// Row step: Ax_new = A x; z = prox_{sigma g*}(z + sigma (2 Ax_new - Ax)); accumulate averages. par = {tau, sigma}.
template <int G, class T>
__global__ void k_row(int m, const int* ptr, const int* idx, const T* val, const T* x, T* Ax, T* z, const T* rl,
                      const T* ru, const T* par, T* Axs, T* zs) {
    int t = blockIdx.x * blockDim.x + threadIdx.x, i = t / G, lane = t % G;
    bool ok = i < m;
    T s = group_dot<G>(ok, i, lane, ptr, idx, val, x);
    if (ok && lane == 0) {
        T sig = par[1], v = z[i] + sig * (2 * s - Ax[i]), q = v / sig;
        T zn = q < rl[i] ? v - sig * rl[i] : (q > ru[i] ? v - sig * ru[i] : T(0));
        Ax[i] = s, z[i] = zn, Axs[i] += s, zs[i] += zn;
    }
}

// Column step: ATz = A^T z, accumulate; if upd, also the next x = proj(x - tau (c + ATz)) for the same column.
template <int G, class T>
__global__ void k_col(int n, const int* ptr, const int* idx, const T* val, const T* z, T* ATz, T* ATzs, int upd,
                      T* x, const T* c, const T* l, const T* u, const T* par, T* xs) {
    int t = blockIdx.x * blockDim.x + threadIdx.x, j = t / G, lane = t % G;
    bool ok = j < n;
    T s = group_dot<G>(ok, j, lane, ptr, idx, val, z);
    if (ok && lane == 0) {
        ATz[j] = s, ATzs[j] += s;
        if (upd) {
            T xn = clampv(x[j] - par[0] * (c[j] + s), l[j], u[j]);
            x[j] = xn, xs[j] += xn;
        }
    }
}

template <class T>
__global__ void k_xupd(int n, T* x, const T* ATz, const T* c, const T* l, const T* u, const T* par, T* xs) {
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) {
        T xn = clampv(x[j] - par[0] * (c[j] + ATz[j]), l[j], u[j]);
        x[j] = xn, xs[j] += xn;
    }
}

// Block reduction of K doubles into part[blockIdx * K + k].
template <int K>
__device__ void block_store(double (&v)[K], double* part) {
    __shared__ double sh[K][32];
    int lane = threadIdx.x % 32, w = threadIdx.x / 32;
    for (int k = 0; k < K; ++k) {
        for (int off = 16; off > 0; off /= 2) v[k] += __shfl_down_sync(0xffffffffu, v[k], off);
        if (lane == 0) sh[k][w] = v[k];
    }
    __syncthreads();
    if (threadIdx.x < 32) {
        int nw = blockDim.x / 32;
        for (int k = 0; k < K; ++k) {
            double s = lane < nw ? sh[k][lane] : 0;
            for (int off = 16; off > 0; off /= 2) s += __shfl_down_sync(0xffffffffu, s, off);
            if (lane == 0) part[blockIdx.x * K + k] = s;
        }
    }
}

template <class T>
__global__ void k_eval_row(int m, double f, const T* Ax, const T* z, const T* rl, const T* ru, const T* Dr,
                           const T* zlast, double* part) {
    double v[NROW] = {};
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < m; i += gridDim.x * blockDim.x) {
        double ax = Ax[i] * f, zz = z[i] * f, lo = rl[i], up = ru[i];
        double r = ax < lo ? lo - ax : (ax > up ? ax - up : 0), ro = r / double(Dr[i]);
        double rr = ro * row_rel_weight(lo, up, Dr[i]);
        v[PR2S] += r * r, v[PR2O] += ro * ro, v[PR2R] += rr * rr;
        if (zz > 0 && up < INFINITY) v[DROW] -= zz * up;
        if (zz < 0 && lo > -INFINITY) v[DROW] -= zz * lo;
        double dz = zz - zlast[i];
        v[DZ2] += dz * dz;
    }
    block_store<NROW>(v, part);
}

template <class T>
__global__ void k_eval_col(int n, double f, const T* x, const T* ATz, const T* c, const T* l, const T* u,
                           const T* Dc, const T* xlast, double* part) {
    double v[NCOL] = {};
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) {
        double xx = x[j] * f, cj = c[j], d = cj + ATz[j] * f, lo = l[j], up = u[j], res = 0;
        v[POBJ] += cj * xx;
        if (d > 0) {
            if (lo > -INFINITY) v[DCOL] += d * lo;
            else res = d;
        } else if (d < 0) {
            if (up < INFINITY) v[DCOL] += d * up;
            else res = -d;
        }
        double ro = res / double(Dc[j]), dx = xx - xlast[j];
        v[DR2S] += res * res, v[DR2O] += ro * ro, v[DX2] += dx * dx;
    }
    block_store<NCOL>(v, part);
}

template <class T>
__global__ void k_sumsq(int n, const T* x, double* part) {  // power iteration norm
    double v[1] = {0};
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) v[0] += double(x[j]) * x[j];
    block_store<1>(v, part);
}

template <class T>
__global__ void k_scale(int n, T* x, double f) {
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) x[j] = T(x[j] * f);
}

__global__ void k_finish(int blocks, int K, const double* part, double* out) {  // fixed-order final sum
    int k = threadIdx.x;
    if (k >= K) return;
    double s = 0;
    for (int b = 0; b < blocks; ++b) s += part[b * K + k];
    out[k] = s;
}

template <class T>
__global__ void k_restart(int n, int avg, double f, T* x, T* xs, T* xlast, T* ys) {  // x <- avg or current
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) {
        T v = avg ? T(xs[j] * f) : x[j];
        x[j] = v, xlast[j] = v, xs[j] = 0, ys[j] = 0;
    }
}

template <class T>
__global__ void k_scale_out(int n, double f, const T* src, double* dst) {
    for (int j = blockIdx.x * blockDim.x + threadIdx.x; j < n; j += gridDim.x * blockDim.x) dst[j] = double(src[j]) * f;
}

static int group_for(const Csr& a) {  // threads per row: the power of two >= mean nnz per row, in [2, 32]
    double mean = a.rows ? double(a.val.size()) / a.rows : 1;
    int g = 2;
    while (g < 32 && g < mean) g *= 2;
    return g;
}

template <class T>
struct Gpu {
    struct Mat {
        int rows = 0, g = 32;
        int *ptr = nullptr, *idx = nullptr;
        T* val = nullptr;
    };
    int m = 0, n = 0, eval_freq = 64;
    bool use_graph = true;
    Mat A, AT;
    T *c, *l, *u, *rl, *ru, *Dr, *Dc, *par;
    T *x, *z, *Ax, *ATz, *xs, *zs, *Axs, *ATzs, *xlast, *zlast;
    T *pv, *pw, *pu;    // power-iteration vectors
    double *part, *out, *hbuf;  // reduction buffers; download buffer
    int red_blocks = 1;
    cudaStream_t st;
    cudaGraphExec_t graph = nullptr;
    std::vector<void*> allocs;
    size_t mem_used = 0;
    std::string name;

    template <class U>
    U* alloc(size_t count) {
        U* p;
        CHECK(cudaMalloc(&p, std::max<size_t>(count, 1) * sizeof(U)));
        allocs.push_back(p);
        return p;
    }
    T* upload(const std::vector<double>& v) {
        std::vector<T> h(v.begin(), v.end());
        T* p = alloc<T>(h.size());
        CHECK(cudaMemcpy(p, h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice));
        return p;
    }
    Mat upload(const Csr& a) {
        Mat d;
        d.rows = a.rows, d.g = group_for(a);
        d.ptr = alloc<int>(a.ptr.size()), d.idx = alloc<int>(a.idx.size());
        CHECK(cudaMemcpy(d.ptr, a.ptr.data(), a.ptr.size() * sizeof(int), cudaMemcpyHostToDevice));
        CHECK(cudaMemcpy(d.idx, a.idx.data(), a.idx.size() * sizeof(int), cudaMemcpyHostToDevice));
        d.val = upload(a.val);
        return d;
    }
    T* zeros(int k) {
        T* p = alloc<T>(k);
        CHECK(cudaMemset(p, 0, std::max(k, 1) * sizeof(T)));
        return p;
    }

    void init(const Scaled& P, int freq, bool graph_on) {
        m = P.m, n = P.n, eval_freq = freq, use_graph = graph_on;
        size_t free0, free1, total;
        CHECK(cudaMemGetInfo(&free0, &total));
        cudaDeviceProp prop;
        CHECK(cudaGetDeviceProperties(&prop, 0));
        name = prop.name;
        CHECK(cudaStreamCreateWithFlags(&st, cudaStreamNonBlocking));
        A = upload(P.A), AT = upload(P.AT);
        c = upload(P.c), l = upload(P.l), u = upload(P.u), rl = upload(P.rl), ru = upload(P.ru);
        Dr = upload(P.Dr), Dc = upload(P.Dc);
        par = zeros(2);
        x = zeros(n), z = zeros(m), Ax = zeros(m), ATz = zeros(n), xs = zeros(n), zs = zeros(m), Axs = zeros(m),
        ATzs = zeros(n), xlast = zeros(n), zlast = zeros(m);
        red_blocks = int(std::min<long>(512, std::max<long>(1, (std::max(m, n) + 255) / 256)));
        part = alloc<double>(size_t(red_blocks) * (NROW + NCOL)), out = alloc<double>(NROW + NCOL);
        pv = zeros(n), pw = zeros(m), pu = zeros(n), hbuf = alloc<double>(std::max(m, n));
        CHECK(cudaDeviceSynchronize());
        CHECK(cudaMemGetInfo(&free1, &total));
        mem_used = free0 - free1;
    }
    ~Gpu() {
        if (graph) cudaGraphExecDestroy(graph);
        for (void* p : allocs) cudaFree(p);
    }

    int eblocks(int k) const { return int(std::min<long>(1024, std::max<long>(1, (k + 255) / 256))); }

    template <class F>
    void dispatch(int g, F&& f) {  // call f with the compile-time group size
        switch (g) {
            case 2: f(std::integral_constant<int, 2>()); break;
            case 4: f(std::integral_constant<int, 4>()); break;
            case 8: f(std::integral_constant<int, 8>()); break;
            case 16: f(std::integral_constant<int, 16>()); break;
            default: f(std::integral_constant<int, 32>());
        }
    }
    int blocks_for(const Mat& a) const { return int((size_t(a.rows) * a.g + 255) / 256); }

    void spmv(const Mat& a, const T* in, T* outv) {
        dispatch(a.g, [&](auto G) {
            k_spmv<decltype(G)::value, T><<<std::max(1, blocks_for(a)), 256, 0, st>>>(a.rows, a.ptr, a.idx, a.val, in, outv);
        });
    }
    double sumsq(int k, const T* v) {
        int b = std::min(red_blocks, eblocks(k));
        k_sumsq<T><<<b, 256, 0, st>>>(k, v, part);
        k_finish<<<1, 32, 0, st>>>(b, 1, part, out);
        double h;
        CHECK(cudaMemcpyAsync(&h, out, sizeof(double), cudaMemcpyDeviceToHost, st));
        CHECK(cudaStreamSynchronize(st));
        return h;
    }
    // Power iteration on A_s^T A_s: est = sqrt(||A^T A v||) with ||v|| = 1 converges to sigma_max(A_s).
    double power_norm(int iters) {
        std::vector<double> h(n);
        double nn = 0;
        for (int j = 0; j < n; ++j) h[j] = start_vec(j), nn += h[j] * h[j];
        std::vector<T> ht(n);
        for (int j = 0; j < n; ++j) ht[j] = T(h[j] / std::sqrt(nn));
        h2d(pv, ht.data(), n * sizeof(T));
        double est = 0;
        for (int it = 0; it < iters; ++it) {
            spmv(A, pv, pw);
            spmv(AT, pw, pu);
            double s = sumsq(n, pu), prev = est;
            if (s <= 0) return 0;
            est = std::sqrt(std::sqrt(s));
            k_scale<T><<<eblocks(n), 256, 0, st>>>(n, pu, 1 / std::sqrt(s));
            std::swap(pv, pu);
            if (it >= 10 && std::fabs(est - prev) <= 1e-7 * est) break;
        }
        return est;
    }

    void h2d(void* dst, const void* src, size_t bytes) {  // ordered on the solver stream (non-blocking stream)
        CHECK(cudaMemcpyAsync(dst, src, bytes, cudaMemcpyHostToDevice, st));
        CHECK(cudaStreamSynchronize(st));
    }
    void set_steps(double tau, double sigma) {
        T h[2] = {T(tau), T(sigma)};
        h2d(par, h, sizeof(h));
    }
    void products() {  // Ax = A x, ATz = A^T z, exactly
        spmv(A, x, Ax);
        spmv(AT, z, ATz);
    }
    void start(const Scaled& P) {
        std::vector<double> x0(n);
        for (int j = 0; j < n; ++j) x0[j] = std::min(std::max(0.0, P.l[j]), P.u[j]);
        std::vector<T> h(x0.begin(), x0.end());
        h2d(x, h.data(), n * sizeof(T));
        h2d(xlast, h.data(), n * sizeof(T));
        products();
        CHECK(cudaStreamSynchronize(st));
    }
    void launch_block(int steps) {
        k_xupd<<<eblocks(n), 256, 0, st>>>(n, x, ATz, c, l, u, par, xs);
        for (int s = 0; s < steps; ++s) {
            dispatch(A.g, [&](auto G) {
                k_row<decltype(G)::value, T><<<std::max(1, blocks_for(A)), 256, 0, st>>>(m, A.ptr, A.idx, A.val, x, Ax, z, rl, ru, par, Axs, zs);
            });
            dispatch(AT.g, [&](auto G) {
                k_col<decltype(G)::value, T><<<std::max(1, blocks_for(AT)), 256, 0, st>>>(
                    n, AT.ptr, AT.idx, AT.val, z, ATz, ATzs, s + 1 < steps, x, c, l, u, par, xs);
            });
        }
    }
    void iterate(int steps) {
        if (use_graph && steps == eval_freq) {
            if (!graph) {  // capture one block of eval_freq iterations once; steps sizes live in device memory
                cudaGraph_t g;
                CHECK(cudaStreamBeginCapture(st, cudaStreamCaptureModeThreadLocal));
                launch_block(steps);
                CHECK(cudaStreamEndCapture(st, &g));
                CHECK(cudaGraphInstantiate(&graph, g, 0));
                CHECK(cudaGraphDestroy(g));
            }
            CHECK(cudaGraphLaunch(graph, st));
        } else {
            launch_block(steps);
        }
        CHECK(cudaGetLastError());
    }
    Kkt eval(bool avg, double f) {
        int br = std::min(red_blocks, eblocks(m)), bc = std::min(red_blocks, eblocks(n));
        k_eval_row<T><<<br, 256, 0, st>>>(m, f, avg ? Axs : Ax, avg ? zs : z, rl, ru, Dr, zlast, part);
        k_finish<<<1, 32, 0, st>>>(br, NROW, part, out);
        k_eval_col<T><<<bc, 256, 0, st>>>(n, f, avg ? xs : x, avg ? ATzs : ATz, c, l, u, Dc, xlast, part + size_t(br) * NROW);
        k_finish<<<1, 32, 0, st>>>(bc, NCOL, part + size_t(br) * NROW, out + NROW);
        double h[NROW + NCOL];
        CHECK(cudaMemcpyAsync(h, out, sizeof(h), cudaMemcpyDeviceToHost, st));
        CHECK(cudaStreamSynchronize(st));
        Kkt k;
        std::copy(h, h + NROW, k.r), std::copy(h + NROW, h + NROW + NCOL, k.c);
        return k;
    }
    void restart(bool avg, double f) {
        k_restart<T><<<eblocks(n), 256, 0, st>>>(n, avg, f, x, xs, xlast, ATzs);
        k_restart<T><<<eblocks(m), 256, 0, st>>>(m, avg, f, z, zs, zlast, Axs);
        products();
    }
    void download(bool avg, double f, std::vector<double>& hx, std::vector<double>& hz) {
        hx.resize(n), hz.resize(m);
        k_scale_out<T><<<eblocks(n), 256, 0, st>>>(n, avg ? f : 1.0, avg ? xs : x, hbuf);
        CHECK(cudaMemcpyAsync(hx.data(), hbuf, n * sizeof(double), cudaMemcpyDeviceToHost, st));
        CHECK(cudaStreamSynchronize(st));
        k_scale_out<T><<<eblocks(m), 256, 0, st>>>(m, avg ? f : 1.0, avg ? zs : z, hbuf);
        CHECK(cudaMemcpyAsync(hz.data(), hbuf, m * sizeof(double), cudaMemcpyDeviceToHost, st));
        CHECK(cudaStreamSynchronize(st));
    }
};


// ---------------------------------------------------------------- CPU backend (same math, std::thread pool)
static inline void cpu_relax() {
#if defined(__x86_64__)
    __builtin_ia32_pause();
#endif
}

// Persistent workers; run(f) calls f(t) for t in [0, nt) and returns when all are done. Workers spin briefly, then
// block on a condition variable, so an oversubscribed machine (other jobs running) degrades instead of collapsing.
class Pool {
public:
    explicit Pool(int nt) : nt_(nt) {
        for (int t = 1; t < nt_; ++t) workers_.emplace_back([this, t] { loop(t); });
    }
    ~Pool() {
        quit_ = true;
        wake();
        for (auto& w : workers_) w.join();
    }
    int size() const { return nt_; }
    void run(const std::function<void(int)>& f) {
        if (nt_ == 1) return f(0);
        job_ = &f;
        done_.store(0);
        wake();
        f(0);
        for (int spins = 0; done_.load() < nt_ - 1; ++spins)
            if (spins > 2000) std::this_thread::yield();
            else cpu_relax();
    }

private:
    void wake() {
        gen_.fetch_add(1);
        if (sleepers_.load() > 0) {
            std::lock_guard<std::mutex> lk(mu_);
            cv_.notify_all();
        }
    }
    void loop(int t) {
        unsigned seen = 0;
        for (;;) {
            int spins = 0;
            while (gen_.load() == seen && ++spins < 2000) cpu_relax();
            if (gen_.load() == seen) {
                std::unique_lock<std::mutex> lk(mu_);
                sleepers_.fetch_add(1);
                cv_.wait(lk, [&] { return gen_.load() != seen; });
                sleepers_.fetch_sub(1);
            }
            seen = gen_.load();
            if (quit_) return;
            (*job_)(t);
            done_.fetch_add(1);
        }
    }
    int nt_;
    std::vector<std::thread> workers_;
    std::atomic<unsigned> gen_{0};
    std::atomic<int> done_{0}, sleepers_{0};
    std::atomic<bool> quit_{false};
    std::mutex mu_;
    std::condition_variable cv_;
    const std::function<void(int)>* job_ = nullptr;
};

template <class T>
struct Cpu {
    struct Mat {
        std::vector<int> ptr, idx, part;  // part: thread boundaries balanced on rows + nnz
        std::vector<T> val;
    };
    int m = 0, n = 0;
    Pool pool;
    Mat A, AT;
    std::vector<T> c, l, u, rl, ru, Dr, Dc, x, z, Ax, ATz, xs, zs, Axs, ATzs, xlast, zlast;
    T tau = 0, sigma = 0;
    std::string name = "cpu";
    size_t mem_used = 0;

    explicit Cpu(int threads) : pool(threads) {}
    static std::vector<T> cvt(const std::vector<double>& v) { return std::vector<T>(v.begin(), v.end()); }
    Mat load(const Csr& a, int nt) {
        Mat d;
        d.ptr = a.ptr, d.idx = a.idx, d.val = cvt(a.val);
        long total = long(a.val.size()) + a.rows;
        d.part.assign(nt + 1, a.rows);
        d.part[0] = 0;
        for (int t = 1; t < nt; ++t) {  // first row whose prefix (rows + nnz) reaches t/nt of the total
            long target = total * t / nt;
            int lo = 0, hi = a.rows;
            while (lo < hi) {
                int mid = (lo + hi) / 2;
                if (long(a.ptr[mid]) + mid < target) lo = mid + 1;
                else hi = mid;
            }
            d.part[t] = lo;
        }
        return d;
    }
    void init(const Scaled& P, int, bool) {
        m = P.m, n = P.n;
        A = load(P.A, pool.size()), AT = load(P.AT, pool.size());
        c = cvt(P.c), l = cvt(P.l), u = cvt(P.u), rl = cvt(P.rl), ru = cvt(P.ru), Dr = cvt(P.Dr), Dc = cvt(P.Dc);
        for (auto* v : {&x, &ATz, &xs, &ATzs, &xlast}) v->assign(n, 0);
        for (auto* v : {&z, &Ax, &zs, &Axs, &zlast}) v->assign(m, 0);
        mem_used = (A.val.size() + AT.val.size()) * (sizeof(T) + sizeof(int)) + size_t(12 * n + 10 * m) * sizeof(T);
    }
    static T dot(const Mat& a, int r, const T* v) {
        T s = 0;
        for (int k = a.ptr[r]; k < a.ptr[r + 1]; ++k) s += a.val[k] * v[a.idx[k]];
        return s;
    }
    void spmv(const Mat& a, const T* in, T* out) {
        pool.run([&](int t) {
            for (int r = a.part[t]; r < a.part[t + 1]; ++r) out[r] = dot(a, r, in);
        });
    }
    double power_norm(int iters) {
        std::vector<double> h(n);
        double nn = 0;
        for (int j = 0; j < n; ++j) h[j] = start_vec(j), nn += h[j] * h[j];
        std::vector<T> pv(n), pw(m), pu(n);
        for (int j = 0; j < n; ++j) pv[j] = T(h[j] / std::sqrt(nn));
        double est = 0;
        for (int it = 0; it < iters; ++it) {
            spmv(A, pv.data(), pw.data());
            spmv(AT, pw.data(), pu.data());
            double s = 0, prev = est;
            for (T v : pu) s += double(v) * v;
            if (s <= 0) return 0;
            est = std::sqrt(std::sqrt(s));
            for (T& v : pu) v = T(v * (1 / std::sqrt(s)));
            std::swap(pv, pu);
            if (it >= 10 && std::fabs(est - prev) <= 1e-7 * est) break;
        }
        return est;
    }
    void set_steps(double t, double s) { tau = T(t), sigma = T(s); }
    void products() {
        spmv(A, x.data(), Ax.data());
        spmv(AT, z.data(), ATz.data());
    }
    void start(const Scaled& P) {
        for (int j = 0; j < n; ++j) x[j] = xlast[j] = T(std::min(std::max(0.0, P.l[j]), P.u[j]));
        products();
    }
    static T clampv(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }
    void xupd(int j) {
        T xn = clampv(x[j] - tau * (c[j] + ATz[j]), l[j], u[j]);
        x[j] = xn, xs[j] += xn;
    }
    void iterate(int steps) {
        // Column phase of step s fused with the x update of step s+1 (same columns, same thread): 2 barriers/step.
        pool.run([&](int t) {
            for (int j = AT.part[t]; j < AT.part[t + 1]; ++j) xupd(j);
        });
        for (int s = 0; s < steps; ++s) {
            pool.run([&](int t) {
                for (int i = A.part[t]; i < A.part[t + 1]; ++i) {
                    T ax = dot(A, i, x.data()), v = z[i] + sigma * (2 * ax - Ax[i]), q = v / sigma;
                    T zn = q < rl[i] ? v - sigma * rl[i] : (q > ru[i] ? v - sigma * ru[i] : T(0));
                    Ax[i] = ax, z[i] = zn, Axs[i] += ax, zs[i] += zn;
                }
            });
            bool upd = s + 1 < steps;
            pool.run([&](int t) {
                for (int j = AT.part[t]; j < AT.part[t + 1]; ++j) {
                    T v = dot(AT, j, z.data());
                    ATz[j] = v, ATzs[j] += v;
                    if (upd) xupd(j);
                }
            });
        }
    }
    Kkt eval(bool avg, double f) {
        int nt = pool.size();
        std::vector<double> pr(size_t(nt) * NROW, 0), pc(size_t(nt) * NCOL, 0);
        const T *Axv = avg ? Axs.data() : Ax.data(), *zv = avg ? zs.data() : z.data();
        const T *xv = avg ? xs.data() : x.data(), *ATzv = avg ? ATzs.data() : ATz.data();
        pool.run([&](int t) {
            double* v = &pr[size_t(t) * NROW];
            for (int i = A.part[t]; i < A.part[t + 1]; ++i) {
                double ax = Axv[i] * f, zz = zv[i] * f, lo = rl[i], up = ru[i];
                double r = ax < lo ? lo - ax : (ax > up ? ax - up : 0), ro = r / double(Dr[i]);
                double rr = ro * row_rel_weight(lo, up, Dr[i]);
                v[PR2S] += r * r, v[PR2O] += ro * ro, v[PR2R] += rr * rr;
                if (zz > 0 && up < kInf) v[DROW] -= zz * up;
                if (zz < 0 && lo > -kInf) v[DROW] -= zz * lo;
                double dz = zz - zlast[i];
                v[DZ2] += dz * dz;
            }
            double* w = &pc[size_t(t) * NCOL];
            for (int j = AT.part[t]; j < AT.part[t + 1]; ++j) {
                double xx = xv[j] * f, cj = c[j], d = cj + ATzv[j] * f, lo = l[j], up = u[j], res = 0;
                w[POBJ] += cj * xx;
                if (d > 0) {
                    if (lo > -kInf) w[DCOL] += d * lo;
                    else res = d;
                } else if (d < 0) {
                    if (up < kInf) w[DCOL] += d * up;
                    else res = -d;
                }
                double ro = res / double(Dc[j]), dx = xx - xlast[j];
                w[DR2S] += res * res, w[DR2O] += ro * ro, w[DX2] += dx * dx;
            }
        });
        Kkt k;
        for (int t = 0; t < nt; ++t) {
            for (int q = 0; q < NROW; ++q) k.r[q] += pr[size_t(t) * NROW + q];
            for (int q = 0; q < NCOL; ++q) k.c[q] += pc[size_t(t) * NCOL + q];
        }
        return k;
    }
    void restart(bool avg, double f) {
        for (int j = 0; j < n; ++j) x[j] = xlast[j] = avg ? T(xs[j] * f) : x[j], xs[j] = 0, ATzs[j] = 0;
        for (int i = 0; i < m; ++i) z[i] = zlast[i] = avg ? T(zs[i] * f) : z[i], zs[i] = 0, Axs[i] = 0;
        products();
    }
    void download(bool avg, double f, std::vector<double>& hx, std::vector<double>& hz) {
        hx.resize(n), hz.resize(m);
        for (int j = 0; j < n; ++j) hx[j] = avg ? double(xs[j]) * f : double(x[j]);
        for (int i = 0; i < m; ++i) hz[i] = avg ? double(zs[i]) * f : double(z[i]);
    }
};

// ---------------------------------------------------------------- the algorithm, shared by both backends
struct Options {
    double tol = 1e-4, time_limit = 3600;
    long max_iter = 2000000000L;
    int eval_freq = 64;
    double row_rel = 0;  // opt-in (0 = off, the old stop test only; 100 for pds-100/fome21): also require ||(row violation)/(1+|row bound|)||_2 <= row_rel * tol 
    bool verbose = false, graph = true;
};

struct Outcome {
    std::string status = "iteration_limit";
    long iterations = 0;
    int restarts = 0;
    bool from_average = false;
    double rel_primal = 0, rel_dual = 0, rel_gap = 0, norm_a = 0, omega = 1, power_s = 0, iterate_s = 0;
    std::vector<double> x, z;  // scaled-space solution
};

template <class B>
Outcome run(B& b, const Scaled& P, const Options& o, double time_budget) {
    Outcome r;
    auto t0 = Clock::now();
    r.norm_a = b.power_norm(200);
    r.power_s = since(t0);
    double eta = r.norm_a > 0 ? 0.95 / r.norm_a : 1.0;
    double cs = 0, bs = 0;
    for (double v : P.c) cs += v * v;
    for (int i = 0; i < P.m; ++i) {
        double lo = P.rl[i], up = P.ru[i];
        double v = std::isfinite(lo) && std::isfinite(up) ? std::max(std::fabs(lo), std::fabs(up))
                   : std::isfinite(lo)                    ? std::fabs(lo)
                   : std::isfinite(up)                    ? std::fabs(up)
                                                          : 0;
        bs += v * v;
    }
    cs = std::sqrt(cs), bs = std::sqrt(bs);
    double w = cs > 1e-10 && bs > 1e-10 ? cs / bs : 1.0;
    b.set_steps(eta / w, eta * w);
    b.start(P);

    auto wkkt = [&](const Kkt& k) {
        double g = k.pobj() - k.dobj();
        return std::sqrt(w * k.r[PR2S] + k.c[DR2S] / w + g * g);
    };
    struct Rel {
        double p, d, g;
        double worst() const { return std::max(p, std::max(d, g)); }
    };
    auto rel = [&](const Kkt& k) {
        double po = k.pobj(), du = k.dobj();
        double prow = o.row_rel > 0 ? std::sqrt(k.r[PR2R]) / o.row_rel : 0;  // per-row relative: no single row can hide in ||b||
        return Rel{std::max(std::sqrt(k.r[PR2O]) / (1 + P.bnorm), prow), std::sqrt(k.c[DR2O]) / (1 + P.cnorm),
                   std::fabs(po - du) / (1 + std::fabs(po) + std::fabs(du))};
    };
    double k_last = wkkt(b.eval(false, 1.0)), k_prev = kInf;
    long it = 0, n_since = 0;
    bool done = false, avg_final = false;
    Rel best{kInf, kInf, kInf};
    while (!done) {
        int steps = int(std::min<long>(o.eval_freq, o.max_iter - it));
        if (steps <= 0) break;
        b.iterate(steps);
        it += steps, n_since += steps;
        double f = 1.0 / double(n_since);
        Kkt kc = b.eval(false, 1.0), ka = b.eval(true, f);
        Rel rc = rel(kc), ra = rel(ka);
        bool use_avg = ra.worst() < rc.worst();
        best = use_avg ? ra : rc, avg_final = use_avg;
        if (o.verbose && (it / o.eval_freq) % 50 == 0)
            std::fprintf(stderr, "it %ld  w %.3e  cur %.2e %.2e %.2e  avg %.2e %.2e %.2e  pobj %.10g\n", it, w, rc.p,
                         rc.d, rc.g, ra.p, ra.d, ra.g, kc.pobj());
        if (!std::isfinite(rc.worst()) && !std::isfinite(ra.worst())) {
            r.status = "numerical_failure";
            break;
        }
        if (best.worst() <= o.tol) {
            r.status = "near_optimal";
            break;
        }
        if (since(t0) > time_budget) {
            r.status = "time_limit";
            break;
        }
        if (it >= o.max_iter) break;
        // Adaptive restart (KKT error with primal weight): to the better of current and average.
        double wc = wkkt(kc), wa = wkkt(ka);
        bool cand_avg = wa < wc;
        const Kkt& kk = cand_avg ? ka : kc;
        double kcand = std::min(wa, wc);
        bool restart = n_since >= 0.36 * it || kcand <= 0.2 * k_last || (kcand <= 0.8 * k_last && kcand > k_prev);
        k_prev = kcand;
        if (restart) {
            double dx = std::sqrt(kk.c[DX2]), dz = std::sqrt(kk.r[DZ2]);
            b.restart(cand_avg, f);
            if (dx > 1e-10 && dz > 1e-10) {
                w = std::exp(0.5 * std::log(dz / dx) + 0.5 * std::log(w));
                b.set_steps(eta / w, eta * w);
            }
            k_last = wkkt(kk), k_prev = kInf, n_since = 0, r.restarts++;
        }
    }
    r.iterations = it;
    r.rel_primal = best.p, r.rel_dual = best.d, r.rel_gap = best.g, r.omega = w, r.from_average = avg_final;
    b.download(avg_final, 1.0 / double(std::max<long>(n_since, 1)), r.x, r.z);
    r.iterate_s = since(t0) - r.power_s;
    return r;
}

// ---------------------------------------------------------------- CLI
static const char* g_dual_path = nullptr;  // --dual OUT: row multipliers for the independent dual-bound check

struct Timings {
    double parse = 0, scale = 0, init = 0, transfer = 0;
};

template <class B>
int finish(B& b, const Model& md, const Scaled& P, const Options& o, const Timings& tm, Clock::time_point t_start,
           const char* json, const char* sol, const char* device, int threads, const char* prec) {
    double setup = tm.parse + tm.scale + tm.init + tm.transfer;
    auto ts = Clock::now();
    Outcome r = run(b, P, o, o.time_limit - setup);
    double solve = since(ts);
    // Unscale and check the point on the ORIGINAL model in fp64.
    std::vector<double> x(P.n);
    for (int j = 0; j < P.n; ++j) x[j] = r.x[j] * P.Dc[j];
    long double obj = md.obj_const;
    for (int j = 0; j < P.n; ++j) obj += (long double)md.cost[j] * x[j];
    std::vector<double> ax(P.m, 0);
    double bviol = 0, rviol = 0;
    for (int j = 0; j < P.n; ++j) {
        for (const Entry& e : md.cols[j]) ax[e.index] += e.value * x[j];
        bviol = std::max(bviol, std::max(md.col_lo[j] - x[j], x[j] - md.col_up[j]));
    }
    for (int i = 0; i < P.m; ++i) rviol = std::max(rviol, std::max(md.row_lo[i] - ax[i], ax[i] - md.row_up[i]));
    double total = since(t_start);
    std::printf("status %s objective %.12g iterations %ld restarts %d rel_p %.2e rel_d %.2e rel_gap %.2e setup %.3fs solve %.3fs total %.3fs device %s\n",
                r.status.c_str(), double(obj), r.iterations, r.restarts, r.rel_primal, r.rel_dual, r.rel_gap, setup,
                solve, total, b.name.c_str());
    if (json)
        if (std::FILE* f = std::fopen(json, "w")) {
            std::fprintf(f,
                         "{\"solver\": \"taral-pdhg (prototype, near-optimal)\", \"status\": \"%s\", \"objective\": %.17g, "
                         "\"rel_primal\": %.6e, \"rel_dual\": %.6e, \"rel_gap\": %.6e, \"tol\": %.3g, \"iterations\": %ld, "
                         "\"restarts\": %d, \"point\": \"%s\", \"max_row_violation\": %.6e, \"max_bound_violation\": %.6e, "
                         "\"setup_s\": %.6f, \"parse_s\": %.6f, \"scale_s\": %.6f, \"device_init_s\": %.6f, \"transfer_s\": %.6f, "
                         "\"solve_s\": %.6f, \"power_iter_s\": %.6f, \"total_s\": %.6f, \"device\": \"%s\", \"backend\": \"%s\", "
                         "\"threads\": %d, \"precision\": \"%s\", \"rows\": %d, \"cols\": %d, \"nnz\": %zu, \"norm_a\": %.6e, "
                         "\"primal_weight\": %.6e, \"eval_freq\": %d, \"graph\": %s, \"device_mem_bytes\": %zu}\n",
                         r.status.c_str(), double(obj), r.rel_primal, r.rel_dual, r.rel_gap, o.tol, r.iterations,
                         r.restarts, r.from_average ? "average" : "current", rviol, bviol, setup, tm.parse, tm.scale,
                         tm.init, tm.transfer, solve, r.power_s, total, b.name.c_str(), device, threads, prec, P.m, P.n,
                         P.A.val.size(), r.norm_a, r.omega, o.eval_freq,
                         (o.graph && !std::strcmp(device, "gpu")) ? "true" : "false", b.mem_used);
            std::fclose(f);
        }
    if (sol)
        if (std::FILE* f = std::fopen(sol, "w")) {
            for (int j = 0; j < P.n; ++j) std::fprintf(f, "%s %.17g\n", md.col_names[j].c_str(), x[j]);
            std::fclose(f);
        }
    // Row multipliers in original row scale (y = Dr * z, minimisation form): any y gives a weak-duality bound that
    // the benchmark harness recomputes independently of this program (gpu/vector_check.py).
    if (g_dual_path)
        if (std::FILE* f = std::fopen(g_dual_path, "w")) {
            for (int i = 0; i < P.m; ++i) std::fprintf(f, "%s %.17g\n", md.row_names[i].c_str(), P.Dr[i] * r.z[i]);
            std::fclose(f);
        }
    return 0;
}

template <class T>
int solve_with(const std::string& device, int threads, const Model& md, const Scaled& P, const Options& o, Timings tm,
               Clock::time_point t_start, const char* json, const char* sol) {
    const char* prec = sizeof(T) == 8 ? "fp64" : "fp32";
    if (device == "gpu") {
        auto t = Clock::now();
        CHECK(cudaFree(nullptr));  // context creation, counted in setup
        tm.init = since(t);
        t = Clock::now();
        Gpu<T> b;
        b.init(P, o.eval_freq, o.graph);
        tm.transfer = since(t);
        return finish(b, md, P, o, tm, t_start, json, sol, "gpu", 0, prec);
    }
    auto t = Clock::now();
    Cpu<T> b(threads);
    b.init(P, o.eval_freq, o.graph);
    tm.transfer = since(t);  // CPU: thread start + copy into solver arrays
    char nm[64];
    std::snprintf(nm, sizeof nm, "cpu x%d threads", threads);
    b.name = nm;
    return finish(b, md, P, o, tm, t_start, json, sol, "cpu", threads, prec);
}

int main(int argc, char** argv) {
    const char *model = nullptr, *json = nullptr, *sol = nullptr;
    std::string device = "gpu";
    int threads = 0;
    bool fp32 = false;
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&] {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "missing value for %s\n", a.c_str());
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--device") device = next();
        else if (a == "--threads") threads = std::atoi(next());
        else if (a == "--tol") o.tol = std::atof(next());
        else if (a == "--time-limit") o.time_limit = std::atof(next());
        else if (a == "--max-iter") o.max_iter = std::atol(next());
        else if (a == "--row-rel") o.row_rel = std::atof(next());
        else if (a == "--eval-freq") o.eval_freq = std::max(1, std::atoi(next()));
        else if (a == "--json") json = next();
        else if (a == "--sol") sol = next();
        else if (a == "--dual") g_dual_path = next();
        else if (a == "--fp32") fp32 = true;
        else if (a == "--no-graph") o.graph = false;
        else if (a == "--verbose") o.verbose = true;
        else model = argv[i];
    }
    if (!model || (device != "gpu" && device != "cpu")) {
        std::fprintf(stderr, "usage: pdhg MODEL.mps --device gpu|cpu [--threads N] [--tol 1e-4] [--time-limit S] "
                             "[--max-iter N] [--json OUT.json] [--sol OUT.sol] [--fp32] [--eval-freq 64] [--row-rel 0] [--no-graph]\n");
        return 2;
    }
    auto t_start = Clock::now();
    Timings tm;
    Model md;
    try {
        md = read_mps(model);
    } catch (const ParseError& e) {
        std::printf("status parse_error: %s\n", e.what());
        return 1;
    }
    tm.parse = since(t_start);
    auto t = Clock::now();
    Scaled P = build_scaled(md, md.maximize ? -1.0 : 1.0);
    tm.scale = since(t);
    if (threads <= 0) {  // auto: measured crossover (see out/gpu ledgers); small models run fastest on one thread
        int hw = int(std::max(1u, std::thread::hardware_concurrency()));
        threads = P.A.val.size() < 200000 ? 1 : hw;
    }
    return fp32 ? solve_with<float>(device, threads, md, P, o, tm, t_start, json, sol)
                : solve_with<double>(device, threads, md, P, o, tm, t_start, json, sol);
}
