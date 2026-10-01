// TARAL-LP GPU track G3: own CSR sparse matrix-vector product kernel (CUDA runtime only, no cuSPARSE).
// One warp per row with a shuffle reduction. The benchmark checks every result against a plain C++ CPU loop
// and times GPU compute against the faster CPU baseline of one thread or every hardware thread (std::thread),
// timing GPU compute (cudaEvent, matrix resident on the device) and host-device transfer separately.
// Build: nvcc -O3 -arch=sm_120 -o csr_matvec gpu/csr_matvec.cu      Run: ./csr_matvec > results.json
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>

#define CHECK(call)                                                                       \
    do {                                                                                  \
        cudaError_t e = (call);                                                           \
        if (e != cudaSuccess) {                                                           \
            std::fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); \
            std::exit(1);                                                                 \
        }                                                                                 \
    } while (0)

template <class T>
__global__ void csr_matvec(int n, const int* __restrict__ ptr, const int* __restrict__ idx, const T* __restrict__ val,
                           const T* __restrict__ x, T* __restrict__ y) {
    int row = (blockIdx.x * blockDim.x + threadIdx.x) / 32, lane = threadIdx.x % 32;
    if (row >= n) return;
    T s = 0;
    for (int k = ptr[row] + lane; k < ptr[row + 1]; k += 32) s += val[k] * x[idx[k]];
    for (int off = 16; off > 0; off /= 2) s += __shfl_down_sync(0xffffffffu, s, off);
    if (lane == 0) y[row] = s;
}

template <class T>
void bench(int n, double density, int repeats, bool last) {
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> u(-1, 1);
    std::bernoulli_distribution keep(density);
    std::vector<int> ptr(n + 1, 0), idx;
    std::vector<T> val, x(n), y_cpu(n), y_gpu(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            if (keep(rng)) idx.push_back(j), val.push_back(T(u(rng)));
        ptr[i + 1] = int(idx.size());
    }
    for (T& v : x) v = T(u(rng));
    size_t nnz = idx.size();

    int nt = int(std::max(1u, std::thread::hardware_concurrency()));
    auto rows = [&](int lo, int hi) {
        for (int i = lo; i < hi; ++i) {
            T s = 0;
            for (int k = ptr[i]; k < ptr[i + 1]; ++k) s += val[k] * x[idx[k]];
            y_cpu[i] = s;
        }
    };
    auto threaded = [&] {  // all hardware threads, contiguous row blocks
        std::vector<std::thread> pool;
        for (int t = 0; t < nt; ++t) pool.emplace_back(rows, int(int64_t(n) * t / nt), int(int64_t(n) * (t + 1) / nt));
        for (auto& th : pool) th.join();
    };
    auto median_time = [&](auto&& fn) {
        fn();  // warm-up
        std::vector<double> ts;
        for (int r = 0; r < repeats; ++r) {
            auto t0 = std::chrono::steady_clock::now();
            fn();
            ts.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }
        std::sort(ts.begin(), ts.end());
        return ts[ts.size() / 2];
    };
    // Honest CPU baseline: the faster of one thread and all threads (thread start-up dominates small sizes).
    double c1 = median_time([&] { rows(0, n); }), cn = median_time(threaded);
    int used = c1 <= cn ? 1 : nt;

    int *dp, *di;
    T *dv, *dx, *dy;
    auto t0 = std::chrono::steady_clock::now();
    CHECK(cudaMalloc(&dp, (n + 1) * sizeof(int)));
    CHECK(cudaMalloc(&di, nnz * sizeof(int)));
    CHECK(cudaMalloc(&dv, nnz * sizeof(T)));
    CHECK(cudaMalloc(&dx, n * sizeof(T)));
    CHECK(cudaMalloc(&dy, n * sizeof(T)));
    CHECK(cudaMemcpy(dp, ptr.data(), (n + 1) * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(di, idx.data(), nnz * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dv, val.data(), nnz * sizeof(T), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dx, x.data(), n * sizeof(T), cudaMemcpyHostToDevice));
    CHECK(cudaDeviceSynchronize());
    double transfer = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    int threads = 256, blocks = int((size_t(n) * 32 + threads - 1) / threads);
    csr_matvec<T><<<blocks, threads>>>(n, dp, di, dv, dx, dy);  // warm-up
    CHECK(cudaGetLastError());
    cudaEvent_t a, b;
    CHECK(cudaEventCreate(&a));
    CHECK(cudaEventCreate(&b));
    std::vector<double> gt;
    for (int r = 0; r < repeats; ++r) {
        CHECK(cudaEventRecord(a));
        csr_matvec<T><<<blocks, threads>>>(n, dp, di, dv, dx, dy);
        CHECK(cudaEventRecord(b));
        CHECK(cudaEventSynchronize(b));
        float ms;
        CHECK(cudaEventElapsedTime(&ms, a, b));
        gt.push_back(ms / 1e3);
    }
    CHECK(cudaMemcpy(y_gpu.data(), dy, n * sizeof(T), cudaMemcpyDeviceToHost));
    double err = 0, scale = 0;
    for (int i = 0; i < n; ++i) err = std::max(err, double(std::abs(y_gpu[i] - y_cpu[i]))), scale = std::max(scale, double(std::abs(y_cpu[i])));
    for (void* p : {(void*)dp, (void*)di, (void*)dv, (void*)dx, (void*)dy}) CHECK(cudaFree(p));

    std::sort(gt.begin(), gt.end());
    double cm = std::min(c1, cn), gm = gt[gt.size() / 2];
    std::printf("  {\"n\": %d, \"dtype\": \"%s\", \"nnz\": %zu, \"cpu_median_s\": %.3e, \"gpu_median_s\": %.3e, "
                "\"cpu_threads\": %d, \"transfer_s\": %.3e, \"speedup_compute\": %.3f, \"speedup_end_to_end\": %.3f, \"max_rel_error\": %.3e}%s\n",
                n, sizeof(T) == 8 ? "float64" : "float32", nnz, cm, gm, used, transfer, cm / gm, cm / (gm + transfer),
                err / std::max(scale, 1e-300), last ? "" : ",");
}

int main() {
    cudaDeviceProp p;
    CHECK(cudaGetDeviceProperties(&p, 0));
    std::printf("{\"tool\": \"taral_g3_csr_matvec\", \"kernel\": \"own warp-per-row CSR, no cuSPARSE\", \"gpu\": \"%s\", "
                "\"cc\": \"%d.%d\", \"density\": 0.01, \"repeats\": 20, \"results\": [\n", p.name, p.major, p.minor);
    int sizes[] = {1024, 2048, 4096, 8192, 16384};
    for (int n : sizes) bench<float>(n, 0.01, 20, false);
    for (int n : sizes) bench<double>(n, 0.01, 20, n == 16384);
    std::printf("]}\n");
}
