# T4 industrial LP measurement

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

Raw result labels: PASS is accepted; PASS_REF or pass* is accepted with an independently checked reference discrepancy; NOANS/noans means no verified answer; PARTIAL means a weaker conclusion than the expected result; BORDER means curvature is within rounding uncertainty; FAIL_STATUS/FAIL_OBJ/FAIL_FEAS/FAIL_KKT mean a status/objective/feasibility/optimality-check failure; REF_UNRESOLVED means no usable reference answer. JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


Single-run CUDA PDHG measurements, fp64, shared Colab free Tesla T4 (15 GB), compute capability 7.5; driver 580.82.07, nvcc 13.0.88, two-vCPU host. Engine cap 300s. Reported source: main the evaluated source snapshot `58f77af` CPU snapshot plus the then-current `gpu/pdhg.cu`; exact GPU source hash was not included with this CSV.

`total_s` includes parsing. `near_optimal` is approximate, not a certificate of an exact optimum. Recheck uses highspy as a reader on the original expanded MPS. HiGHS objective references are absent for both rail cases. A tolerance in the filename describes the solver stopping rule, not a bound on objective error or every original-model row residual.

At 1e-4, PDS-100 objective relative error is 7.72e-4; FOME21 and PDS-100 original-model relative row violations are 1.23 and 1.02. At 1e-6 they are still 0.0108 and 0.0309, not Netlib-strict feasibility. Both rail cases reach the 300s cap at 1e-6. These measurements do not establish a speed advantage over exact simplex.

[Byte-exact received CSV](ledger.csv).
