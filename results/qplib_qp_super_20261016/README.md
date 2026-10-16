# QPLIB_8602 run record - opt-in supernodal QP kernel (October 16, 2026)

This folder records the run behind the README claim that QPLIB_8602 reached engine-reported optimal with the opt-in `TARAL_QP_SUPER=1` supernodal LDL^T kernel, raising the [historical QPLIB selection](../qplib_widening_20261003/) tally to 13/19 reaching engine-reported optimal. Engine-reported optimal with independently verified primal feasibility is not a rigorous optimality certificate.

## What was run

- Model: QPLIB_8602 from the QPLIB library zip (qplib.zib.de), download sha256 `cb452aeece377ea346095498267ed30fb184a432820f980a22441b603db45231`; expanded to MPS (sha256 `ef338fd9b14a9fd84d178242d609653eb3458fb75fca27dacbee181e43d919b6`) with the conversion script used in the audit.
- Engine: checked-in `main` at merge `c1363ef` (October 16), built twice from clean trees with g++ 11.4.0, `g++ -O2 -std=c++17 -pthread src/*.cpp`; both builds byte-identical, sha256 `616e709e82d0945f175da1074dfe66033cfbd8caac78bf6df97f7be892a571c0`.
- Flags: `--method ipm --time-limit 300`, environment `TARAL_QP_SUPER=1` (default 2e8-flop cutoff).

## Result

Status `optimal`, 14 iterations, objective **27353.387205546067**, dual objective 27353.387192656512, primal residual 1.009e-12, dual residual 0.0, relative gap 4.712e-10, wall 112.5 s on the audit host (shared sandbox; the engine's own host reported 101 s). Raw JSON: [8602_qp_super.json](8602_qp_super.json); returned point: [8602_qp_super.sol](8602_qp_super.sol).

## Independent checks

- **Scalar path, same checked-in source, knob unset** at a 900-second cap: status `optimal`, 14 iterations, objective 27353.387205546063, wall 474.1 s - the same objective to 4e-12 absolute. Raw JSON: [8602_scalar_900s.json](8602_scalar_900s.json). The supernodal path is roughly 4x faster on wall here.
- **Independent KKT check** of the returned point against the original `.qplib` file ([kkt_check_8602.txt](kkt_check_8602.txt)): objective confirmed 27353.38720554609; max row violation 1.009e-12, bound violation 0; Q symmetric (symmetry error 0). The checker's from-scratch LSQR dual estimate did not converge tightly (relative residual 1.0e-2, istop 2), so the dual side stands on the engine's self-reported dual residual 0.0 and gap 4.712e-10; the primal side is independently verified.
- **HiGHS cross-check**: a highspy solve of the same expanded MPS was attempted on the audit host and had not completed after 40 minutes, so no HiGHS reference objective is included; the independent checks above stand without it.
- **17-case identity sample** (GENHS28, HS21, HS35, HS51, CVXQP1_S, DUAL1, DUALC1, AUG3DQP, CONT-050, GOULDQP3, QSCFXM3, ic0/7/13/21/33/39): status, objective and iterations identical three ways (pre-merge build, merged build knob unset, merged build `TARAL_QP_SUPER=1`), 17/17.
- **Fallback exercised**: on QPLIB_10038 with `TARAL_QP_SUPER=1` at the default cutoff the supernodal factor goes non-finite; the engine logs `ipm: supernodal factor non-finite, scalar fallback` when `TARAL_QP_SUPER_VERBOSE` is set and continues on the scalar kernel, reaching the same `time_limit` trajectory as the unset path (3 iterations, identical primal residual at a 60-second cap) instead of the `numerical_failure` the pre-fallback build produced. A quick regression for the kernel lives at [tests/qp_super_regression.sh](../../tests/qp_super_regression.sh).

## Ledger row

[ledger.csv](ledger.csv) uses the historical selection's column format. QPLIB_8500 remains unsolved (time limit), consistent with the historical ledger.
