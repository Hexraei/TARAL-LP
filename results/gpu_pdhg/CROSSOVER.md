# PDHG-to-simplex crossover prototype (local; prototype, not a headline claim)

## Round 2 (measured Oct 2, after the LU speed-up and the factorization deadline; this section supersedes the numbers below)
Same protocol as round 1: `gpu/pdhg.cu` at tol 1e-6 writes the primal point and row multipliers, `taral --method dual --warm-sol X.sol
--warm-dual X.dual --cross-tol 1e-5` builds the basis (src/crossover.cpp, unchanged) and the dual simplex finishes exactly.
Cold dual = same binary, same session, no flags. Medians of 3 alternating runs on a quiet machine (RTX 5050 Laptop + Intel Core 7 240H).
Total = PDHG process time + simplex wall time.

| Case | PDHG (solve / process) | Warm basis factor | Simplex iterations warm / cold | Simplex wall warm / cold | Total warm vs cold |
|---|---|---|---|---|---|
| dfl001 (6,071 rows) | 0.59 s / 2.16 s (0.68 s on the one run with a warm card) | first 0.12 s; 2.98 s in 251 refactorizations over the run (cold: 1.62 s in 179) | 29,869 / 21,256 | 18.20 s / 11.88 s | 20.4 s vs **11.88 s** (best case 18.9 s) |
| dfl001, cross-tol 1e-3 | same | 0.12 s first | 29,474 / 21,256 | 17.60 s / 11.88 s | 19.8 s vs 11.88 s |
| packing_100000 (50,000 rows) | 1.74 s / 2.23 s | **not factored in 29.9 s** (see below) | 30,128 (after fallback, in the remaining time) / 35,160 | no answer in 120 s / no answer in 120 s | no answer either way; the cold dual objective is not reached |

All finished dfl001 routes reached the same objective as the cold dual (11266396.0467). Cold dual dfl001 runs: 11.75, 11.88, 12.05 s.

Factorization deadline (task 1 of this round): `SparseLU::factor` now takes an optional deadline. A timeout surfaces as `time_limit` (exit 4), never as a
singular basis. A warm basis gets a quarter of the limit to factor; if it does not, the solver logs the fallback (stderr and the JSON/status message,
e.g. `warm basis not factored within 29.875480 s (after 45175 of 50000 pivots, 21291882 entries left in the active submatrix), slack basis used`)
and continues from the slack basis. The previously hanging packing_100000 handoff now returns at the limit (60.07 s with a 60 s limit, 120.01 s
with 120 s; 1.6 GB resident at the stop) instead of running past 11.5 minutes.

### Why the warm start loses
dfl001: the starting basis built from the PDHG point is poor on both sides. Of 6,071 basis slots only 3,507 are determined by variables farther than 1e-5 from a bound
(3,474 at PDHG 1e-8; the PDHG point sits inside the optimal face, it is not a vertex). The basis the rest of the slots produce has 2,585 dual infeasibilities (so the
dual simplex has to run its phase 1: 12,083 of the 29,869 iterations) and 4,421 primal infeasible basic variables (sum of infeasibilities 3.0e5). Phase 2
then needs 17.8k iterations against 21.3k from scratch, a 16% saving that the 12k phase 1 iterations more than eat. The PDHG point is accurate enough (rel_p 1e-6, rel_d 4e-10, objective
error 5e-7): the failure is basis selection, not point quality. Tolerance sweep with this build (single runs, cross-tol 1e-5):

| PDHG tol | PDHG process | variables off bound | start: dual infeasibilities / primal infeasible basics | phase 1 its | total its | simplex wall |
|---|---|---|---|---|---|---|
| 1e-4 | 1.93 s | 3,573 | 2,434 / 4,413 | 11,805 | 23,653 | 12.7 s |
| 1e-6 | 2.16 s | 3,507 | 2,585 / 4,421 | 12,083 | 29,869 | 18.2 s |
| 1e-8 | 2.28 s | 3,474 | 2,503 / 4,417 | 13,282 | 30,532 | 17.0 s |

A tighter tolerance does not improve the basis; none of the three beats the cold dual (11.88 s) even before PDHG time is added. One basis-selection variant was tried:
ranking variables by the ratio of distance-to-bound to normalised reduced cost (a complementarity score) instead of distance first, reduced cost for leftovers. It was much worse
(3,455 dual infeasibilities, primal infeasibility sum 3.6e6, no finish in 60 s) and was not kept.

packing_100000: this is not a tolerance problem. The vertex-like basis (50,140 variables off a bound for 50,000 slots, same at PDHG 1e-8) has a **dense nucleus**: when
the factorization is stopped, 45,175 of 50,000 pivots are done and 21.3M entries remain in the active submatrix, i.e. the last ~4,800 pivots work on a
nearly dense ~4,800 x 4,800 block (4,800^2 = 23M). The random sparse structure (5 entries per column in random rows) makes any basis of this shape fill in. At PDHG 1e-8 (19.4 s process time) the picture
is the same (45,347 pivots, 20.7M entries at the stop). Even if the factorization finished (a dense switch for the nucleus is not implemented; its cost is not measured), every iteration would carry L and U of
order 10^7 entries: an estimate (not measured) of 20 to 40 ms per ftran/btran pair, against 3.4 ms per iteration for the cold dual (35,160 iterations in 120 s), so the warm start
could win only by saving most of the iterations, which dfl001 gives no reason to expect.

### Round 2 conclusion
The handoff still does not pay on either case: dfl001 20.4 s against 11.9 s cold, packing_100000 no exact answer either way. The LU speed-up made the cold dual
faster (14.1 s to 11.1-11.9 s), which narrows the room a handoff could use. PDHG alone remains the fast route to a near-optimal answer (dfl001 5e-7 objective error in 0.6 s of solve time;
packing_100000 in 1.7 s of solve time). A useful crossover needs a basis-identification step (push/pull from the PDHG point to a vertex), not a better ranking of the
same point, and for the packing family a factorization that tolerates a dense nucleus.

## Round 1 (measured earlier on the previous LU, kept for the record)

Handoff: `gpu/pdhg.cu` writes an approximate primal point and row multipliers; `taral --warm-sol X.sol --warm-dual X.dual
--cross-tol T` (src/crossover.cpp) turns them into a simplex basis and the existing warm-start path finishes exactly.
Basis rule: the m variables (structurals and row logicals) farthest from their bounds are basic; slots left over when fewer
than m are farther than T from a bound go to the smallest reduced costs |c - A'y|; everything else sits at its nearest bound.
Hardware: RTX 5050 Laptop (PDHG), Intel Core 7 240H (simplex). Default engine path is unchanged without the flags.

## dfl001 (6,071 rows; the exact primal simplex does not finish in 60 s)
| Route | PDHG | Simplex after handoff | Total | Simplex iterations |
|---|---|---|---|---|
| Cold dual simplex (no handoff) | - | 11.8 s | **11.8 s** | 21,256 |
| PDHG tol 1e-4 + dual simplex (cross-tol 1e-5) | 0.47 s (1.9 s on a first run) | 12.7 s | 13.2 s | 23,653 |
| PDHG tol 1e-6 + dual simplex | 0.68 s | 17.8 s | 18.5 s | 29,869 |
| PDHG tol 1e-8 + dual simplex | 3.7 s | 18.7 s | 22.4 s | 30,532 |
| PDHG tol 1e-6 + primal simplex | 0.68 s | not done in 60 s | - | - |
| PDHG tol 1e-6, primal point only (no multipliers), dual or primal | 0.68 s | not done in 60 s | - | - |
All handoff routes that finished reached the same objective as the cold dual, 11266396.0467 (HiGHS agrees).

Finding: the handoff does **not** pay on dfl001. The PDHG point is not near a vertex: only about 3,500 of the 6,071 basis slots
are determined by variables clearly away from a bound, and the rest is guesswork from reduced costs, so the simplex still
needs about as many iterations as from the slack basis (23.7k against 21.3k) and the PDHG time is added on top. A tighter
PDHG tolerance makes the basis worse here, not better. (The cold dual used to take 31.8 s on the Kaggle T4 host and 14-15 s
here before the LU speed-up; it is 11.8 s here now, which is why the PDHG advantage of 0.5 s approximate versus 12 s exact
does not translate into an exact-answer speedup.) PDHG alone remains the fast route to a near-optimal dfl001
(objective error 4e-5 at 1e-4, 5e-7 at 1e-6 against HiGHS), separate protocol line.

## packing_100000 (50,000 rows, 100,000 columns, 500k nonzeros; synthetic, gpu/gen_lp.py)
- PDHG tol 1e-6: 2.3 s, objective -68881.0220 against HiGHS (interior point) -68881.0219 in 2,331 s.
- Cold dual simplex: no answer in 120 s (34,656 iterations).
- Handoff: the warm basis was still inside its **initial LU factorization after 11.5 minutes (1.9 GB resident)** and was killed.
  The 50,000-row basis built from a PDHG point has a dense nucleus. The factorization has no time check, so the 120 s limit
  was not honoured: a known defect of the prototype path (a basis that cannot be factored in reasonable time should fall back
  to the slack basis; not done).
Conclusion: no crossover speed-up demonstrated on either case. This is the honest outcome; a usable crossover needs a
real basis-identification step (push/pull to a vertex), not a nearest-to-bound crash.
