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

## a) Correctness (600 required cases: seeds 26119 3 7 26120 31 47, 100 each)

Commit 0 turns the three `parse_error` cases (seed 26120 RND63 = 87, RND64 = -165; seed 31 RND94 = 163) into matches.
Every ledger below: **600/600 match HiGHS, 0 wrong, 0 reference-side cases** (no HiGHS-infeasible/taral-feasible,
no taral-better-than-reference); 565 optimal, 35 infeasible; status and objective identical to `ac46975` on all 600.
Two of the 600 have no integer column (seed 26119 RND65, seed 3 RND49): they go through the LP path and report no node count.

| binary | total nodes | wrong |
|---|---|---|
| `ac46975` (before) | 3344 | 0 |
| `31aa19f` propagation, A (default) | 2713 | 0 |
| `31aa19f` `--no-prop-prune` (B) | 2722 | 0 |
| `2346d62` + reduced-cost fixing, A (default) | **2423** | 0 |
| `2346d62` `--no-prop-prune` (B) | 2432 | 0 |

Per seed, before → after (`2346d62` A): 26119 448→307, 3 630→492, 7 548→425, 26120 559→381, 31 639→437, 47 520→381.
Instances with fewer nodes: 247; more: 6 (largest RND27 of seed 26119, 6 → 13; none changes a result).
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
The six LP-feasible cases (one node each): seed 26119 RND90, seed 3 RND98, seed 26120 RND68, seed 31 RND1 and RND52, seed 47 RND50.
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
python benchmarks/milp_ledger.py random --engine ./taral --commit HASH --out results/milp_propagation/random_HASH_A.csv      # seeds 26119 3 7 26120 31 47, 100 each, 60 s
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
