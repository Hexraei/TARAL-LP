# Measured published LP applications

Evaluated source snapshot `34ea8b4`; solver sources clean at measurement time.
Warning-free GCC build, 30-second limit, one local run per case. Source and
binary hashes, build flags, model hashes, complete variable points and
independent checks are in `results.json`. Model translation and original
citations are in [the example guide](../../examples/literature_lp/).

| Application | Engine result | Independent reference | Original row/bound violation |
| --- | --- | --- | --- |
| Steel production | optimal, 192000 illustrative dollars/week, 3 iterations | HiGHS 192000 | 0 / 0 |
| Transportation | optimal, 153.675 thousand dollars, 7 iterations | HiGHS 153.675 | 0 / 0 |

Relative objective discrepancy is zero in both stored comparisons. The checker
re-parses each original MPS independently of the C++ reader and recomputes
activity and objective from the returned point. HiGHS supplies a separate
optimal reference. This is a floating-point comparison, not exact proof.
Millisecond process times are retained for completeness, not speed claims.
No industrial outcome or MILP application coverage is established here.
