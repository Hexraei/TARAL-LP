# MIPLIB widening: 40 measured cases

## Reading the technical terms

LP means linear programming; MILP means mixed-integer linear programming; QP means quadratic programming. IPM is an interior-point method. PDHG is primal-dual hybrid gradient, the approximate GPU method. KKT (Karush-Kuhn-Tucker) checks test feasibility and optimality conditions. A "gate" is the stated validation rule, not an exact-arithmetic proof. "Strict" means the separately stated tighter tolerance; a non-strict pass meets the ordinary rule but not that tighter check. A ledger is a per-case result table. Source hashes, file paths and command flags are retained only so engineers can reproduce a measurement. Historical measurements are not current-source claims.

JSON is the machine-readable result format. Warm starts reuse a previous solution or simplex basis; cold starts do not. fp64 is double-precision floating point; FMA is fused multiply-add.


The raw `ledger.csv` records six objective matches, seven infeasibility matches,
26 engine time limits and one no-JSON run (`neos-1425699`). Zero rows are labeled
wrong under this protocol. That is not a claim of 40 successful solves or proof
that the run without a machine-readable result is harmless.

Evaluated source snapshot `58f77af`, not the later Q integration, was reported for this Kaggle run.
The reported notebook is `hexraei/notebook89d59928c7`, version 1. Environment,
build/source provenance and CSV hash are in `provenance.json`; per-instance
download and MPS hashes are in the CSV. The MPS bytes and per-case logs were
not delivered with this table, so no byte-level replay or crash adjudication is
claimed here.

The engine requested 120 seconds and HiGHS 300 seconds, soft caps (maximum
recorded walls 120.53 and 300.04 seconds). The current implementation times out
on 26 measured cases. HiGHS also reached its own cap on 11 of those: `ej`,
`p2m2p1m1p0n100`, `markshare_5_0`, `markshare1`, `gen-ip054`, `gen-ip016`,
`neos-5140963-mincio`, `gen-ip002`, `gen-ip021`, `k16x240b`,
`neos-3046615-murg`. These unequal caps do not establish a same-budget speed
comparison.

`neos-1425699` returned no JSON after 50.53 seconds; the reference reported an
optimal objective 3179698977 in 0.42 seconds. The cause remains unknown. A crash
is a hypothesis under investigation, not a result proved by the CSV.
