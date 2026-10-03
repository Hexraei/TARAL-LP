# MILP bound propagation and reduced-cost fixing (branch `cloud/milp-propagation`)

Base `442ca16`. Cloud container (Intel Xeon 2.8 GHz, 4 cores). **Cloud timings are indicative only; the
claim is the node reduction, not speed. Nothing here goes into a data pack without a Kaggle re-run.**
Build line unchanged and warning-free at every commit:
`g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp`.
Not touched: `src/lu.cpp dual.cpp simplex.cpp ipm.cpp`, the README, notebooks, the Netlib corpus/tolerances/denominator.

## Commits (one change each, gated before the next)

| commit | change |
|---|---|
| `ac46975` | 0. MPS reader: BOUNDS naming a column absent from COLUMNS creates a zero-cost, entry-less column (HiGHS accepts it) |
| `9a66928` | tools: `benchmarks/milp_ledger.py`, `benchmarks/milp_stress.py` |
| `31aa19f` | 1. activity-based bound propagation on integer columns (root and after each branching) |
| `2346d62` | 2. reduced-cost fixing at nodes with an incumbent |

Reference stack here: scipy 1.17.1 / highspy 1.15.1 (the brief mentions scipy 1.15.3; its seed-31 RND39 false-infeasible did not appear).

## a) Correctness (600 required cases: seeds 26000+119 3 7 26000+120 31 47, 100 each)

Commit 0 turns the three `parse_error` cases (seed 26000+120 RND63 = 87, RND64 = -165; seed 31 RND94 = 163) into matches.
Every ledger below: **600/600 match HiGHS, 0 wrong, 0 reference-side cases** (no HiGHS-infeasible/taral-feasible,
no taral-better-than-reference); 565 optimal, 35 infeasible; status and objective identical to `ac46975` on all 600.
Two of the 600 have no integer column (seed 26000+119 RND65, seed 3 RND49): they go through the LP path and report no node count.

| binary | total nodes | wrong |
|---|---|---|
| `ac46975` (before) | 3344 | 0 |
| `31aa19f` propagation, A (default) | 2713 | 0 |
| `31aa19f` `--no-prop-prune` (B) | 2722 | 0 |
| `2346d62` + reduced-cost fixing, A (default) | **2423** | 0 |
| `2346d62` `--no-prop-prune` (B) | 2432 | 0 |

Per seed, before → after (`2346d62` A): 26000+119 448→307, 3 630→492, 7 548→425, 26000+120 559→381, 31 639→437, 47 520→381.
Instances with fewer nodes: 247; more: 6 (largest RND27 of seed 26000+119, 6 → 13; none changes a result).
Supplementary, not an acceptance set: 2000 further cases (seeds 1001..1020 × 100), 0 wrong, nodes 11080 → 9005 → 8062;
1000 larger models with equality rows (`benchmarks/milp_stress.py --seeds 1..10`), 0 wrong, nodes 208620 → 187959 → 159814.

### Pruning rule and its audit (program-owner ruling, Option A default)

Propagation prunes a node only when an **integer** column's domain is empty, `ceil(lb - eps) > floor(ub + eps)`
with `eps >= 1e-6` (plus a term growing with the row magnitude, so borderline cases are never pruned). Continuous
columns are never tightened, and a row with no integer column is left to the LP. A node can have a feasible LP
relaxation and still hold no integer point (an integer column with bounds [0.3, 0.7]), so "the LP finds the pruned
node feasible" is not a failure. The audit (`--audit-prop`) instead searches the node's pre-propagation box with a
plain branch and bound (no propagation, no fixing; zero objective, or the incumbent as cutoff for boxes reduced-cost
fixing emptied). It also re-solves the LP, and, for commit 2, searches the region each reduced-cost fixing removes
for any point better than the incumbent. `--no-prop-prune` (Option B) discards the propagated bounds and lets the LP decide.

| set, final binary `2346d62` | pruned by propagation | LP infeasible | LP feasible (integer-empty by plain B&B) | integer-feasible (a bug) | fixings audited | fixings wrong |
|---|---|---|---|---|---|---|
| 600 required | 157 | 150 | 7 | **0** | 1510 | **0** |
| 2000 extra | 479 | 457 | 22 | **0** | 4798 | **0** |
| 1000 stress | 4099 | 3487 | 612 | **0** | 121691 | **0** |

At `31aa19f` (propagation only) on the 600: 56 prunes = 50 LP-infeasible + 6 LP-feasible/integer-empty, 0 integer-feasible.
The six LP-feasible cases (one node each): seed 26000+119 RND90, seed 3 RND98, seed 26000+120 RND68, seed 31 RND1 and RND52, seed 47 RND50.
Reproduce: `python benchmarks/milp_ledger.py random --engine WRAPPER --commit HASH --out OUT.csv` where WRAPPER is a script
running `taral --audit-prop "$@"`; the audit vector in the CSV is `[nodes, lp_infeasible, lp_feasible, lp_other, int_empty,
int_feasible, undecided, cutoff_only, rc_checked, rc_bad]`. The 157 count above is for the final rule; an earlier draft that also pruned on a row-activity
violation reported 163 and had no audit, and is not what ships. `rc_skipped` (basis not trusted) was 0 in every set.

## b) Existing LP gates (93-case Netlib, 60 s, strict; `benchmarks/netlib_gate.py`, same container, HiGHS reference cached from the repo)

| run | `442ca16` | `2346d62` |
|---|---|---|
| A (default `taral`) | 90/93, strict 89/90, 0 wrong | 90/93, strict 89/90, 0 wrong |
| B (`taral --method dual`) | 91/93, strict 90/91, 0 wrong | 91/93, strict 90/91, 0 wrong |

All 94 cases (93 + TRUSS) have identical verdict and strict flag before and after; `pilot.we`/`pilot4` stay the fixed
exclusion, `dfl001` is the Run A time-limit case. Ledgers in `netlib/`. This run overlapped with other jobs, so its
timings are not comparable; the pass/strict counts are.

## c) Benefit on the fixed instance lists (60 s per instance, same lists before and after)

Random: see above. Small-MIPLIB list = the 26 pinned MIPLIB 3 instances of `benchmarks/milp_miplib.py` (the "MIPLIB-easy"
ledger was not in the repo; this pinned list is the closest). These runs were made one at a time with nothing else running
(an earlier overlapped "before" run was discarded because its timings were inflated). Cells are `status objective / nodes / seconds`
(`opt` = optimal, `TL` = time limit; `-` = no incumbent).

| instance | HiGHS | `ac46975` | `31aa19f` A | `31aa19f` B | `2346d62` A |
|---|---|---|---|---|---|
| p0033 | 3089 | opt 3089 / 2,013 / 0.1 | opt 3089 / 551 / 0.0 | opt 3089 / 724 / 0.0 | opt 3089 / 411 / 0.0 |
| p0201 | 7615 | opt 7615 / 1,279 / 2.2 | opt 7615 / 1,158 / 2.0 | opt 7615 / 1,158 / 2.4 | opt 7615 / 985 / 2.0 |
| p0282 | 258411 | opt 258411 / 4,899 / 4.8 | opt 258411 / 3,947 / 3.7 | opt 258411 / 3,947 / 4.3 | opt 258411 / 1,950 / 2.1 |
| p0548 | 8691 | TL - / 17,667 / 60.0 | TL - / 32,138 / 60.0 | TL - / 24,103 / 60.0 | TL - / 32,317 / 60.0 |
| egout | 568.1007 | opt 568.1007 / 68,482 / 22.4 | opt 568.1007 / 68,482 / 21.1 | opt 568.1007 / 68,482 / 22.9 | opt 568.1007 / 68,482 / 23.5 |
| flugpl | 1201500 | opt 1201500 / 7,705 / 0.4 | opt 1201500 / 2,356 / 0.1 | opt 1201500 / 2,870 / 0.2 | opt 1201500 / 2,348 / 0.1 |
| lseu | 1120 | opt 1120 / 33,798 / 6.0 | opt 1120 / 19,260 / 2.7 | opt 1120 / 19,712 / 3.3 | opt 1120 / 8,279 / 1.3 |
| mod008 | 307 | opt 307 / 18,735 / 2.2 | opt 307 / 18,735 / 2.4 | opt 307 / 18,735 / 2.8 | opt 307 / 6,687 / 0.9 |
| stein27 | 18 | opt 18 / 9,797 / 8.0 | opt 18 / 9,458 / 7.0 | opt 18 / 9,458 / 7.3 | opt 18 / 9,484 / 9.0 |
| stein45 | 30 | TL 30 / 17,741 / 60.0 | TL 30 / 19,514 / 60.0 | TL 30 / 18,259 / 60.0 | TL 30 / 16,886 / 60.0 |
| bell5 | 8966406.492 | TL - / 80,593 / 60.1 | TL 8966406.492 / 105,915 / 60.1 | TL 8966406.492 / 107,546 / 60.1 | TL 8966406.492 / 100,913 / 60.1 |
| bell3a | 878430.316 | TL 878430.316 / 60,599 / 60.0 | TL 878430.316 / 60,845 / 60.0 | TL 878430.316 / 58,630 / 60.0 | opt 878430.316 / 38,204 / 37.3 |
| misc03 | 3360 | opt 3360 / 1,740 / 5.0 | opt 3360 / 1,423 / 1.8 | opt 3360 / 1,423 / 2.0 | opt 3360 / 1,183 / 1.8 |
| misc07 | 2810 | TL 2810 / 6,641 / 60.0 | TL 2810 / 14,000 / 60.0 | TL 2810 / 12,240 / 60.0 | TL 2865 / 9,245 / 60.0 |
| enigma | 0 | opt 0 / 913 / 0.2 | opt 0 / 486 / 0.1 | opt 0 / 1,087 / 0.2 | opt 0 / 486 / 0.1 |
| gt2 | 21166 | TL 22342 / 206,996 / 60.2 | opt 21166 / 8,387 / 2.0 | TL 22342 / 175,179 / 60.2 | opt 21166 / 1,644 / 0.6 |
| khb05250 | 106940226 | opt 106940226 / 2,336 / 1.9 | opt 106940226 / 2,336 / 1.9 | opt 106940226 / 2,336 / 2.0 | opt 106940226 / 2,300 / 2.4 |
| vpm1 | 20 | TL 22 / 28,236 / 60.0 | TL 22 / 30,300 / 60.0 | TL 22 / 30,300 / 60.0 | TL 21 / 28,954 / 60.0 |
| vpm2 | 13.75 | TL - / 22,489 / 60.0 | TL - / 24,531 / 60.0 | TL - / 23,085 / 60.0 | TL - / 24,554 / 60.0 |
| dcmulti | 188182 | opt 188182 / 4,079 / 9.5 | opt 188182 / 3,682 / 8.5 | opt 188182 / 3,682 / 8.7 | opt 188182 / 3,792 / 10.7 |
| rgn | 82.2 | opt 82.2 / 5,049 / 1.1 | opt 82.2 / 5,405 / 1.0 | opt 82.2 / 5,405 / 1.0 | opt 82.2 / 3,329 / 0.7 |
| pk1 | 11 | TL 11 / 67,648 / 60.0 | TL 11 / 74,597 / 60.0 | TL 11 / 70,502 / 60.0 | TL 11 / 80,014 / 60.1 |
| pp08a | 7350 | TL 8730 / 44,007 / 60.0 | TL 8730 / 46,737 / 60.0 | TL 8730 / 43,999 / 60.0 | TL 8730 / 43,534 / 60.1 |
| noswot | -41 | TL - / 29,252 / 60.0 | TL - / 32,672 / 60.0 | TL - / 31,081 / 60.0 | TL - / 32,490 / 60.0 |
| fixnet6 | 3983 | TL - / 20,209 / 60.0 | TL - / 19,915 / 60.0 | TL - / 19,523 / 60.0 | TL - / 22,346 / 60.0 |
| blend2 | 7.598985 | opt 7.598985 / 15,284 / 51.6 | opt 7.598985 / 11,878 / 34.8 | opt 7.598985 / 11,878 / 34.5 | opt 7.598985 / 10,267 / 34.0 |

| list | `ac46975` | `31aa19f` A | `31aa19f` B | `2346d62` A |
|---|---|---|---|---|
| solved (optimal) | 14 | 15 | 14 | **16** |
| wrong | 0 | 0 | 0 | 0 |
| nodes, instances solved by both `ac46975` and the column | 176,109 | 149,157 (14) | 150,897 (14) | 119,983 (14) |
| nodes, all 26 (includes time-limited runs) | 778,187 | 618,708 | 765,344 | 551,084 |

Acceptance (c): no instance changes status or objective wrongly (no verdict `WRONG`; every time-limited row keeps bound ≤ HiGHS ≤ incumbent),
and solved count and node totals improve for both commits, so neither is reverted.

### Losses and caveats (shown, not hidden)

- **misc07** ends at the time limit with a worse incumbent at `2346d62` (2865 vs 2810). Still within the honest-limit rule, but it is a regression on that instance.
- Wall time rose on solved instances despite fewer or equal nodes: dcmulti 9.5 → 10.7 s, stein27 8.0 → 9.0 s, khb05250 1.9 → 2.4 s, egout 22.4 → 23.5 s (identical nodes). Reduced-cost fixing refactorises the basis at every node with an incumbent, so per-node cost goes up when fixing finds little. Nothing here supports a speed claim.
- pk1 and pp08a/noswot/fixnet6/p0548/vpm2/stein45 stay unsolved. Their node counts at the limit depend on node speed and are not a measure of progress.
- `31aa19f` Option B does not solve gt2 (Option A does): the direct prune matters there, 175,179 vs 8,387 nodes.
- Propagation's pruning is not LP-verified; the audit above is the evidence, not a proof. Pruning on LP status and the "never prune on Unbounded/NumericalFailure/TimeLimit" rule are unchanged.
- The random models are tiny, so node counts there say little about large models; the MIPLIB list is only 26 instances at 60 s on one container.

## d) Ledgers and exact command lines

Every CSV has a `.meta.json` next to it with the commit, the engine flags, the exact command lines, seeds and tool versions.
`build`: `g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp`

```
python benchmarks/milp_ledger.py random --engine ./taral --commit HASH --out results/milp_propagation/random_HASH_A.csv      # seeds 26000+119 3 7 26000+120 31 47, 100 each, 60 s
python benchmarks/milp_ledger.py random ... --seeds 1001 ... 1020 --ref-cache REF.json                                         # extra2000_*
python benchmarks/milp_ledger.py miplib --engine ./taral --commit HASH --out results/milp_propagation/miplib_HASH_A.csv        # 26 pinned instances, 60 s
python benchmarks/milp_ledger.py compare BEFORE.csv AFTER.csv                                                                   # per-instance losses/gains
python benchmarks/netlib_gate.py --out OUT --engine ./taral --corpus CORPUS --time-limit 60 --with-extra                       # Run A
python benchmarks/netlib_gate.py --out OUT --engine WRAPPER_WITH_--method_dual --corpus CORPUS --time-limit 60 --with-extra    # Run B
```

| ledger | what |
|---|---|
| `random_before_ac46975`, `random_31aa19f_{A,B_no-prop-prune,A_audit}`, `random_2346d62_{A,B_no-prop-prune,A_audit}` | 600 required cases |
| `extra2000_{ac46975,31aa19f_A_audit,2346d62_A_audit}` | supplementary 2000 cases |
| `miplib_{before_ac46975,31aa19f_A,31aa19f_B_no-prop-prune,2346d62_A}` | 26 pinned MIPLIB 3 instances; the 31aa19f/2346d62 runs used binaries built before the commits were created and verified byte-identical (`cmp`) to clean builds of those commits |
| `netlib/ledger_*`, `netlib/summary_*` | 93-case Netlib protocol, `442ca16` and `2346d62`, Run A and Run B |

The stress runs (`benchmarks/milp_stress.py`) print one summary line and are not stored as CSV.

---

# e) Follow-up (C4): the four MIPLIB losses and the big-values `numerical_failure`

Branch tip before this follow-up `05aaaa1`. The hashes in sections a)-d) are the names from before the authorship rewrite:
`ac46975` = `aeec670`, `9a66928` = `f1bc952`, `31aa19f` = `d979350`, `2346d62` = `c3acc2f`. Same container class as before
(4 cores); a faster machine than the one behind the ledgers above, so **node counts at a time limit are not comparable across the
two machines** and every comparison below is made on one machine. The documented 600-case count reproduces exactly here
(`c3acc2f` A = 2,423 nodes). No speed claim is made anywhere in this section.

| commit | change |
|---|---|
| `6bc0ef2` | reduced-cost fixing backs off after calls that fix nothing (skip 1, 3, 7, then up to 15 calls; a fixing or a new incumbent resets it) |
| `61dfed4` | a node whose LP ends neither Optimal nor Infeasible has its widest integer column bisected instead of being left unresolved; adds `cpp-engine/tests/milp/big_values_unresolved.mps` |
| `f2fcfaf` | `--audit-prop` runs the nested plain search on every pruned node, LP-infeasible ones included |
| `9674587` | tool `benchmarks/milp_bigvals.py` (large-magnitude generator) |
| `9c9e7b1` | a node LP that fails on a propagated box is retried on the box without that propagation |
| `4dc330e` | scrub: ledger seeds in the public 26000+N form |

Only `src/milp.cpp` changed in the engine. Tolerance constants (`kIntTol kFeasTol kGapTol kPropEps kRcMin kRcDualTol` and the
LP engine's own checks) are untouched.

## e1) The four documented losses: diagnosis

Everything below is in `followup/miplib_loss_diagnosis.csv`, `followup/reduced_cost_time_share.csv`, `followup/rc_backoff_cap_sweep.csv`.
Builds of `aeec670 d979350 c3acc2f` and the final tree were each compiled three ways (`-O3 -march=native`, the same with
`-ffp-contract=off`, and generic `-O2`); the node counts of a complete run are deterministic for a given build.

**misc07 (worse incumbent at the limit): search-path luck, not a mechanism.** With a 12,000-node limit every build ends at 2810 (the
HiGHS value). The node at which 2810 is first found depends on the floating-point build far more than on the commit:

| commit | native | `-ffp-contract=off` |
|---|---|---|
| `aeec670` (before) | 6,703 | 1,743 |
| `d979350` (propagation) | 11,243 | 5,778 |
| `c3acc2f` (+ reduced-cost fixing) | 3,878 | 9,134 |
| final (`9c9e7b1`) | 3,764 | 9,851 |

Within one commit the spread (up to 6x) exceeds the spread between commits, and the 2865 in the 60 s ledger (9,245 nodes) does
not reproduce here: this container's `c3acc2f` finds 2810 at node 3,878. The incumbent at a time limit on a hard instance is
a lottery over the branching path; nothing in propagation or fixing worsens it systematically. No code change.

**dcmulti (nodes 3,792 vs 4,079 in the ledger, 10.7 s vs 9.5 s): path noise in the nodes, overhead in the time.** Complete-run node
counts (native / `-ffp-contract=off`): before 2,962 / 3,422, propagation 3,566 / 3,390, + fixing 3,227 / 3,517, final 3,088 / 3,398.
A 15% swing between builds of one commit is as large as any commit-to-commit difference, so the node count says nothing here.
The time difference is the reduced-cost refactorization: the timer inside `fix_by_reduced_cost` measured 0.63 s of 5.79 s.

**stein27 (8.0 s -> 9.0 s) and khb05250 (1.9 s -> 2.4 s): per-node overhead of reduced-cost fixing, confirmed.** Nodes are flat
(stein27 9,515 propagation vs 9,551 with fixing; khb05250 2,336 vs 2,300) while the timer measures

| instance | `c3acc2f`: calls, fixings, seconds in the call / total | final: calls, fixings, seconds / total |
|---|---|---|
| stein27 | 4,770, 71, 0.60 / 5.63 | 313, 3, 0.04 / 5.18 |
| khb05250 | 1,215, 1,070, 0.13 / 1.45 | 796, 1,041, 0.08 / 1.41 |
| dcmulti | 1,729, 555, 0.63 / 5.79 | 609, 582, 0.21 / 5.09 |

(a diagnostic timer, other jobs were running; it locates the cost, it is not a benchmark). stein27 fixed a bound in 1.5% of its
calls. This is what `6bc0ef2` changes. It does not change the search for the better in any way I can claim: it removes fixings,
so node counts move a little in both directions on the solved instances (same-machine, `c3acc2f` -> final: p0282 1,995 -> 2,343,
rgn 2,795 -> 3,045, misc03 1,273 -> 819, lseu 8,504 -> 7,991, mod008 6,687 -> 6,367, khb05250 2,300 -> 2,354) and gt2 swings:

gt2 is bimodal. With the back-off cap at 0 (= `c3acc2f`) it solves in 4,583 nodes, with cap 1 in 282,953, with cap 2 and 3 it hits
the limit, with cap 4 (shipped) it solves in 37,752 (`followup/rc_backoff_cap_sweep.csv`); `aeec670` and `d979350` do not solve it at
60 s on this machine. Any perturbation of the branching path flips it, so its `c3acc2f` result was not evidence of a mechanism and
its final result is not evidence against one. Same-machine totals over the 16 instances solved at 60 s
(`followup/miplib_*_solved16_samemachine.csv`): nodes 536,537 (`aeec670`, gt2 unsolved at 301,342), 566,294 (`d979350`),
168,465 (`c3acc2f`), 211,192 (final); 16/16 optimal for the last two. All 26 at 60 s with the final tree
(`followup/miplib_final_A.csv`): 16 optimal, 10 time-limited, **0 wrong**; every time-limited row keeps bound <= HiGHS <= incumbent
(verdict `unsolved` = honest limit; misc07 ends at 2810 in 20,182 nodes here).

## e2) The big-values `numerical_failure`

The reproducer you mentioned is not in the repository or in this session, so I regenerated one of that shape with
`benchmarks/milp_bigvals.py` (row scales up to 1e7, matrix entries up to 9e7, costs 1e3..1e9, integer columns and a known feasible
point). Seed 504, case 17 of the `all` family = `cpp-engine/tests/milp/big_values_unresolved.mps` (15 columns, 12 rows). It
**reaches the true optimum -1,820,000 (HiGHS: -1,820,000) and then ends `numerical_failure`** with one unresolved node: 51 nodes
at `aeec670` (before any propagation), 50 at `d979350`, 46 at `c3acc2f`. So it predates this branch. If your own file is
a different model, send it; it may need another look.

Cause: the node LP message is `final point violates constraints by 0.000000`. `Simplex::finish` (`src/simplex.cpp`) rejects a point
whose scaled violation exceeds 1e-7 and reports NumericalFailure; the value is below 5e-7 (the message prints six digits), within the
MILP layer's 1e-6 check. `solve_lp_dual` and the cold retry fail the same way. The LP engine files are not touched on this branch, so
the fix is in the search: `61dfed4` bisects the widest integer column of such a node. The halves cover the box, each gets an LP, and
nothing is pruned on the failed solve. It ends `optimal -1820000` in 76 nodes (74 with the later `9c9e7b1`), with and without
`--no-prop-prune`, under `--audit-prop`, and agrees with HiGHS. On this model the failure cleared within 5 consecutive bisections;
the cap is 8 (256 leaves), past which, or with no integer column left to split, the node stays unresolved and the run still ends
`numerical_failure` honestly.

Soundness of the fallbacks was fault-injected with scratch builds (not committed): 3% / 10% of node LPs forced to NumericalFailure -
0 wrong optimal/infeasible in either; all 600 still solve at 3%, at 10% four end `numerical_failure` at nodes with no integer column
to split; 30% of propagated nodes forced to fail their first LP (exercises `9c9e7b1`): 600/600 match, audit 0 integer-feasible.

## e3) What the large-magnitude hunt also found (not fixed here)

`milp_bigvals.py` on families `mid wide row all`, seeds 500..515 x 40 (640 models each), final tree, 10 s: 0 wrong
optimal/infeasible claims from the engine's side (the two flagged mismatches are reference-side, see below), 0 `numerical_failure`, 34 time limits (`followup/bigvals_sweep_final.txt`).

* **Propagation exposed LP failures at the root** (`mixed` family, rows scaled 1e-3..1e7): seed 57 case 23 solves in 3 nodes before
  propagation, ended "root LP failed" after `d979350`; `9c9e7b1` solves it again. Seed 57 case 8 does not fail but stalls (next point).
* **LP-engine stall (src/dual.cpp / simplex.cpp, out of scope here).** The 34 time limits: 30 also time out at `aeec670`. Of the four
  that `aeec670` solves in 10 s, `mid` seed 501 case 29 is only slower (153,666 nodes against 100,541; every build solves it in 20 s) and
  three are new stalls (`all` seed 502 case 8 since `d979350`, `row` seed 503 case 37 since `c3acc2f`, `row` seed 511 case 21 in the final tree
  only; `c3acc2f` solves that one in 11,761 nodes). Instrumented: the dual simplex finishes in 2-9 iterations, then the primal clean-up runs 3.6-4.4 million iterations on a
  model of at most 25 columns until the time limit. That is cycling in the LP clean-up, and the path (box) decides whether a run meets it.
  A per-node time cap in `milp.cpp` would be time-dependent and would also kill legitimate long root LPs on big models; not done.
* **Tolerance edge, not changed.** The row tolerance is relative to `1 + |rhs|`, so with right-hand sides of 1e6..1e7 it admits absolute
  violations of a few units. In the `range` family (`--family range`, e.g. seed 168 case 30) the engine reports
  3,262,492 where HiGHS reports 3,262,355.6; the engine's point has relative row violation 6.6e-7 (HiGHS's: 0). Neither is wrong
  under its own tolerance. In two `mid` cases the reference is wrong: seed 501 case 33 (HiGHS "infeasible"; the engine's point is exactly
  feasible and integral) and seed 504 case 39 (the engine's exactly feasible point is better than HiGHS's by 4.6e-5 relative).
* Five `mixed` models fail the root LP identically at `aeec670` (seeds 16, 60, 60, 61, 65): LP engine, not changed.

## e4) Acceptance (final tree = `9c9e7b1` for src/; ledgers in `followup/`)

600 required cases (seeds 26000+119 3 7 26000+120 31 47, 100 each), same reference stack as before:

| binary | total nodes | wrong | status / objective |
|---|---|---|---|
| `c3acc2f` A, this machine | 2,423 | 0 | 565 optimal, 35 infeasible |
| final, A (default) | **2,424** | 0 | identical to `c3acc2f` on all 600 |
| final, B (`--no-prop-prune`) | **2,431** | 0 | identical |

Per seed, A / B: 26000+119 307 / 309, 3 493 / 495, 7 425 / 425, 26000+120 381 / 382, 31 437 / 437, 47 381 / 383. Against `c3acc2f` one
instance differs (seed 3, +1 node); A against B: 5 instances fewer, 1 more. **600/600 match HiGHS, 0 wrong, 0 reference-side.**

Final-tree audit rerun (`--audit-prop`, every propagation-pruned node gets the nested plain search with propagation and fixing
off; vector `[nodes, lp_infeasible, lp_feasible, lp_other, int_empty, int_feasible, undecided, cutoff_only, rc_checked, rc_bad]`):

| set | pruned by propagation | LP infeasible | LP feasible, integer-empty by plain B&B | integer-feasible (a bug) | undecided | fixings audited | fixings wrong |
|---|---|---|---|---|---|---|---|
| 600 required | 154 | 148 | 6 | **0** | 0 | 1,489 | **0** |
| 2000 extra (seeds 1001..1020 x 100; 2000/2000 match) | 482 | 459 | 23 | **0** | 0 | 4,764 | **0** |
| 1000 stress (`milp_stress.py --seeds 1..10`; 0 wrong) | 4,024 | 3,421 | 603 | **0** | 0 | 117,287 | **0** |

The 600 audit run also matches HiGHS on all 600, with the same 2,424 nodes as the plain run (the audit does not alter the search).
The nested search on the LP-infeasible nodes starts from the same LP the audit just found infeasible, so it is a re-check of the LP
engine rather than new evidence; it is run because every pruned node was asked for.

Reproduce: `python benchmarks/milp_ledger.py random --engine WRAPPER --commit HASH --out OUT.csv` (WRAPPER runs `taral --audit-prop "$@"`
or `--no-prop-prune`), `python benchmarks/milp_ledger.py miplib ...`, `python benchmarks/milp_bigvals.py --engine ./taral --family all --seeds 504 --n 18`
(case 17 is the reproducer). Build: `g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp`, warning-free at every commit
(`6bc0ef2 61dfed4 f2fcfaf 9c9e7b1`); each stage binary was byte-compared with the build that was tested.
