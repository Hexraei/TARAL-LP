# Adversarial QP tests and optimality certificates

A seeded, deterministic stress suite for the engine's convex-QP path (`QUADOBJ`/`QMATRIX`/`QSECTION` objective, interior
point), in the same spirit as [adversarial.md](adversarial.md) (LP/MILP). Tools live in `tools/advqp/`, the runner and seed
list in `tests/`. **Tests and tools only: nothing under `src/` was changed**; bugs are reported as reproducers, not fixed.

Engine under test: `src/` at the tip of `cloud/adversarial-tests` (`7d0af4d`), built `g++ -O3 -march=native -std=c++17`.
Reference: HiGHS 1.15.1 (`highspy`) QP solver. Cloud container, 4 cores; the whole suite takes about 3 minutes (engine time
is 8 s of that; the rest is HiGHS and the certificate LPs).

```
tests/run_adversarial_qp.sh               # full suite, 1780 cases
QUICK=1 tests/run_adversarial_qp.sh       # 8 cases per category
SELFTEST=1 tests/run_adversarial_qp.sh    # also the generator/writer/reader self-check (can be slow: HiGHS hangs on a few models)
```

Exit status 0 means every failure is listed in `tests/qp_known_failures.tsv` (id, class, reproducer, detail); 1 means a new
failure or a reproducer regression. Seeds: `tests/qp_seeds.txt` (category *k* uses base `70000 + 1000*k`, case *i* uses
`base + i`; generators use numpy `RandomState` only). Needs `g++`, `python3`, `numpy`, `highspy`.

## What is checked per case

The generator's own numbers are the model (`Q`, `c`, `A`, row/column bounds, sense, constant). The case is written as MPS,
the engine runs with `--sol --json`, and the returned point is judged by code that shares nothing with `src/`:

1. **Exact primal check** (`qpcheck.primal_check`): row and column violations in rational arithmetic on the doubles that
   were written, reported **absolute** and **relative** (`/(1+|violated bound|)`), plus the exact objective
   `c'x + 0.5 x'Qx + const`. The gate is the relative `1e-6`; the absolute value is recorded for every point and
   `abs>1e-6` is counted separately (see the second table).
2. **Optimality certificate from the point alone** (`qpcheck.dual_certificate`), no engine residual used. Two LPs over
   multipliers `y+,y-` (row sides) and `z+,z-` (column bounds), all with the right signs:
   * stationarity: smallest `t` with `|(g - A'(y+ - y-) - (z+ - z-))_j| <= t (1+|g_j|)`, `g = c + Qx` (`stat rel`; the
     unscaled version is `stat abs`);
   * certified gap: among multipliers with that residual (up to rounding), the smallest
     `sum y*(distance of A x to its bound) + z*(distance of x to its bound)`. With exact stationarity this equals
     `f(x) - D(x,y,z)` (Wolfe dual), so for convex `Q` it bounds `f(x) - f*` from above (`gap rel = gap/(1+|f|)`).
   Gate: `stat rel <= 1e-6` and `gap rel <= 1e-6`.
3. **Objective** vs the constructed optimum when the generator knows it, else vs HiGHS (`|obj - target| <= 1e-6 max(1,|target|)`),
   and the engine's *reported* objective must equal the exactly recomputed one.
4. **Truth by construction** where possible, never taken from a solver: most feasible categories build the polyhedron around
   a dyadic point `x0` and choose `c = A'y + z - Q x0` with correctly signed (some zero, i.e. weakly active) multipliers, so
   `x0` is optimal (KKT is sufficient for convex `Q`); infeasible cases have an arithmetic contradiction with margin
   `>= 1e-3` (checked again by an independent phase-1 LP when the engine says infeasible); unbounded cases carry a ray with
   integer data, verified exactly (`qpcheck.verify_ray`: `Ad` in the row cone, `d` in the column cone, `Qd = 0`, `c'd` improves).
   Before a case is judged, the constructed optimum itself must pass the certificate (else `GEN_MISMATCH`, never a pass).
5. **Convexity** is decided from the eigenvalues of `Q` (minimisation sense): clearly PSD (`lam_min >= -1e-13 |Q|`), clearly
   indefinite (`<= -1e-6 |Q|`), otherwise `border`. For clearly indefinite `Q` only `nonconvex` is accepted.
6. **Reference**: HiGHS built from the generator data with `passModel` (not from the MPS file), in a forked child with a hard
   cap; any non-determinate answer is re-solved with presolve off. HiGHS cannot take indefinite `Q`, so none is used there.

Outcome columns: `pass`; `pass*` (engine certified, reference shown wrong); `partial` (`dual_infeasible`: verified recession
ray but primal feasibility not established, correct but weaker than the truth); `border` (`Q` within numerical noise of
semidefinite, any defensible answer with a passing certificate); `wrong`; `noans` (time/iteration limit,
`numerical_failure`, refusal); `other` (reference unresolved, generator or harness error). `rate = (pass + pass*)/cases`;
`strict` additionally requires the exact **absolute** violation `<= 1e-6`.

## Results (full run, 1780 cases)

```
category               cases  pass pass* partial border wrong noans other  known    rate  strict
------------------------------------------------------------------------------------------------
pd_baseline               80    80     0       0      0     0     0     0      0  100.0%  100.0%
psd_rank_def              80    80     0       0      0     0     0     0      0  100.0%  100.0%
psd_free_null             60    56     4       0      0     0     0     0      0  100.0%  100.0%
unbounded                 90    65     0      25      0     0     0     0      0   72.2%   72.2%
infeasible               100   100     0       0      0     0     0     0      0  100.0%  100.0%
indefinite                80    80     0       0      0     0     0     0      0  100.0%  100.0%
indef_convex_feasible     40     0     0       0      0     0    40     0     40    0.0%    0.0%
near_psd                  80    47     4       0     26     0     0     3      3   63.8%   63.8%
zero_q                    60    60     0       0      0     0     0     0      0  100.0%  100.0%
degenerate_cons          100   100     0       0      0     0     0     0      0  100.0%  100.0%
active_set                80    80     0       0      0     0     0     0      0  100.0%  100.0%
bounds_special            90    90     0       0      0     0     0     0      0  100.0%  100.0%
bounds_huge               80    69     0       0      0     1    10     0     11   86.2%   86.2%
bounds_tight              60    60     0       0      0     0     0     0      0  100.0%  100.0%
scale_q                   60    60     0       0      0     0     0     0      0  100.0%  100.0%
scale_cols                60    60     0       0      0     0     0     0      0  100.0%  100.0%
scale_rows                60    59     0       0      0     0     1     0      1   98.3%   96.7%
ill_cond                  60    35     0       0      0     0    25     0     25   58.3%   58.3%
maximize                  60    59     0       0      0     0     1     0      1   98.3%   98.3%
obj_const                 60    60     0       0      0     0     0     0      0  100.0%  100.0%
qform                     60    60     0       0      0     0     0     0      0  100.0%  100.0%
empty_structure           70    70     0       0      0     0     0     0      0  100.0%  100.0%
portfolio                 40    40     0       0      0     0     0     0      0  100.0%  100.0%
large_sparse              20    20     0       0      0     0     0     0      0  100.0%  100.0%
fuzz_mixed               150   150     0       0      0     0     0     0      0  100.0%  100.0%
------------------------------------------------------------------------------------------------
TOTAL                   1780  1640     8      25     26     1    77     3     81   92.6%   92.5%

category               points max abs viol  abs>1e-6  rel>1e-6 max stat rel  max gap rel  max act rel  act>1e-6 max |x-x*|
--------------------------------------------------------------------------------------------------------------------------
pd_baseline                80    9.910e-09         0         0    0.000e+00    9.016e-09    5.790e-01        41   9.77e-04
psd_rank_def               80    1.033e-08         0         0    0.000e+00    9.117e-09    1.679e-02        19   0.00e+00
psd_free_null              60    9.069e-09         0         0    0.000e+00    1.527e-09    0.000e+00         0   0.00e+00
unbounded                   0    0.000e+00         0         0    0.000e+00    0.000e+00    0.000e+00         0   0.00e+00
infeasible                  0    0.000e+00         0         0    0.000e+00    0.000e+00    0.000e+00         0   0.00e+00
indefinite                  0    0.000e+00         0         0    0.000e+00    0.000e+00    0.000e+00         0   0.00e+00
indef_convex_feasible       0    0.000e+00         0         0    0.000e+00    0.000e+00    0.000e+00         0   0.00e+00
near_psd                   61    3.352e-12         0         0    0.000e+00    9.205e-09    4.462e-07         0   0.00e+00
zero_q                     60    6.777e-10         0         0    0.000e+00    2.892e-09    0.000e+00         0   0.00e+00
degenerate_cons           100    1.250e-09         0         0    0.000e+00    9.198e-09    3.061e-03        26   2.71e-04
active_set                 80    4.999e-09         0         0    0.000e+00    9.644e-09    1.740e-03        34   1.04e-03
bounds_special             90    5.270e-10         0         0    1.360e-09    9.494e-09    2.483e-02        30   5.65e-04
bounds_huge                70    1.133e-09         0         0    0.000e+00    1.860e-03    7.183e-01        31   2.50e-02
bounds_tight               60    5.990e-09         0         0    0.000e+00    7.129e-09    7.415e-03        36   6.65e-04
scale_q                    60    7.533e-10         0         0    5.078e-07    8.712e-09    9.991e-01        30   5.97e-02
scale_cols                 60    2.032e-10         0         0    0.000e+00    8.639e-09    5.028e-01        25   6.40e-03
scale_rows                 59    7.688e-06         1         0    0.000e+00    8.118e-09    7.104e-03        24   1.33e-03
ill_cond                   35    3.514e-13         0         0    0.000e+00    3.346e-08    3.920e-05         2   0.00e+00
maximize                   59    1.178e-09         0         0    0.000e+00    9.904e-09    9.910e-03        27   1.06e-03
obj_const                  60    1.405e-08         0         0    0.000e+00    7.464e-09    6.364e-01        29   2.11e-03
qform                      60    6.248e-10         0         0    0.000e+00    7.812e-09    3.399e-03        21   3.23e-04
empty_structure            60    1.368e-11         0         0    0.000e+00    9.577e-09    2.494e-03        10   6.55e-04
portfolio                  40    3.450e-12         0         0    0.000e+00    8.973e-09    2.550e-04         9   0.00e+00
large_sparse               20    6.612e-11         0         0    0.000e+00    6.604e-09    3.297e-03        20   0.00e+00
fuzz_mixed                150    7.242e-09         0         0    8.159e-08    9.661e-09    1.112e-01        55   6.35e-04
```

Reading the first table: `partial` and `border` are not failures but they are not passes either, so the rate understates
neither. 17 of 25 categories are at 100 %. The second table covers only cases where the engine returned a point (1404 of
1780): `abs>1e-6` is the number of points whose exact absolute row/bound violation exceeds `1e-6` (rel gate not hit),
`rel>1e-6` those failing the gate itself, `act>1e-6` the points where the strict *active-set* stationarity (below) exceeds `1e-6`.

### Absolute versus relative violation

Across the 1404 returned points the largest exact absolute violation is `1.0e-8` except in one case:
`scale_rows-86051` has **7.7e-6 absolute** (2.96e-12 relative), a row scaled by about `1e6`, which the relative gate hides.
That is the only point with `abs > 1e-6`; 7 points exceed `1e-8`; 871 points have some nonzero violation. The engine's own
`max_row_viol`/`max_bound_viol` agree with the independent exact absolute value to within `5e-10` on every point (it reports
`7.7e-6` for that case too), so its self-report was honest about absolute violation; only the gate's relative form hides it.

## Failures, honestly

Nothing here was fixed or hidden; every row below has ids in `tests/qp_known_failures.tsv` and a reproducer in `tests/qp_repro/`.

| class | cases | what happens |
| --- | --- | --- |
| `false_optimal_gap` (**wrong**) | 1 | `bounds_huge-82051` (`qp_false_optimal_gap.mps`, 5 columns, 2 rows, a ranged row `[-1e15, 0.625]`): the engine says `optimal`, its own json reports `gap 1.1e-10`, but the exact objective is 33.28062 against the true 33.34375 (relative `1.9e-3`); the independent certificate gives a gap of `0.064`. The only wrong answer in 1780 cases, and the engine's self-reported gap was wrong. Not root-caused. |
| `numfail_huge_bounds` | 10 of 80 `bounds_huge` | `numerical_failure` ("no progress in 30 iterations" or "converged on the scaled problem but not on the original model") on convex QPs whose inactive bounds are `1e6`..`1e19` (finite). Failures occur at every magnitude tried, `1e6` included (2 of 10), and never for bounds `>= 1e20` (infinite by convention, 30 of 30 pass). |
| `numfail_ill_cond` | 25 of 60 `ill_cond` | `numerical_failure` when `cond(Q)` is `1e10` (11 of 15) or `1e12` (14 of 15); all 30 cases at `1e6`/`1e8` pass. |
| `numfail_other` | 2 | `maximize-88059` (n=4, m=5) and `scale_rows-86044` (row scale span 9e3): "no progress in 30 iterations" on small, well-posed models. |
| `refuses_indefinite_q_convex_on_feasible_set` | 40 of 40 | `Q` indefinite but positive definite on the null space of the equality rows (unique optimum, computed from the KKT linear system in numpy): the engine answers `nonconvex`. A deliberate limitation (it needs `Q` PSD), counted as `noans` because the model is convex on its feasible set. |
| `dual_infeasible_on_unbounded` (**partial**) | 25 of 90 `unbounded` | Feasible models with a verified ray spanning two columns (`Ad = 0` through opposing coefficients): the engine reports `dual_infeasible` instead of `unbounded` in 25 of 30 such cases (the unit-ray families all pass, 60 of 60). |
| `reference_error_near_psd` (**other**) | 3 | HiGHS returns an error or "not set" on near-semidefinite models; the engine's point is exactly feasible with zero stationarity residual and gap `<= 1e-8`, but there is no reference objective, so they are not scored as passes. |

### Reference (HiGHS) defects found on the way

* `psd_free_null` (4 cases) and `near_psd` (4 cases), all with no rows (`m = 0`) and free columns: HiGHS returns "Optimal" with objective
  `0` while the gradient is not zero; the engine's objective is lower (e.g. -1184.5) and certified. Scored `pass*`.
* 33 `Solve error`, 27 `not set` and 9 `time limit` statuses from HiGHS; a constructed-unbounded model on which HiGHS returned
  "Optimal" at `x ~ 4.5e7` (found while validating the generator; those cases are judged by construction).
* HiGHS 1.15.1 prints diagnostics to stdout from inside forked solves, which broke the harness protocol once (fixed: the
  worker keeps its JSON on a private descriptor).
* The opt-in self-check found that HiGHS' own MPS reader reads `QSECTION OBJ` listings differently from the engine and from
  this suite's convention (both triangles listed, objective `0.5 x'Qx`); the engine matches the suite on all 60 `qform` cases
  (`QUADOBJ` either token order, `QMATRIX`, `QSECTION OBJ`, split duplicate lines that sum, explicit zero and cancelling lines).
  That is a convention note, not an engine failure.

## A design note: the first certificate was too strict, and why both numbers are kept

The first version of the check let only constraints within `1e-6` of a bound carry multipliers. A smoke run failed almost
every case, including ones whose objective was right to `1e-9`. The reason: an interior-point answer without crossover sits
slightly *inside* bounds that are active at the optimum (for example `x = 0.87533` against a bound `0.875`), so that check
measures distance to the optimal face, not optimality. The gated certificate was therefore changed to the all-multiplier
form above, whose gap is a real bound on `f(x) - f*`. The strict active-set number is kept and **reported, not gated**
(`act rel`, `act>1e-6`): **469 of 1404 points (33 %) exceed `1e-6` and 162 exceed `1e-3`**. Where the optimum is known and
unique (652 points) the distance to it, `max |x - x*|`, has median `3e-7`, exceeds `1e-6` for 298 points and `1e-4` for 194,
worst `6e-2` (flat objectives: tiny curvature). So the engine's objective and duality gap are trustworthy to its stated
`1e-8`, but its coordinates are only about as accurate as the curvature allows; anyone who needs a vertex or the active set
must not read it off the returned `x`. This change was made after seeing the 5-per-category smoke run and before the full run;
the tolerances (`1e-6`) were not changed.

Two harness details a reader should know: the gap LP is given a stationarity allowance equal to the first LP's optimum
(or, when HiGHS reports that infeasible by rounding, the smallest allowance in a half-decade ladder that works), and
distances below `1e-13` are treated as rounding noise. With these, no `GEN_MISMATCH` or `HARNESS_ERR` remains.

## Categories

`pd_baseline`, `psd_rank_def` (rank 0 upward), `psd_free_null` (free columns, minimiser not unique), `unbounded`
(three ray families), `infeasible` (five contradiction families, margins 1e-3..4), `indefinite` (flipped eigenvalue,
negative definite, zero diagonal saddle, small negative direction), `indef_convex_feasible`, `near_psd` (relative
perturbations -1e-4..+1e-10 of a PSD `Q`), `zero_q` (empty/zero/cancelling `QUADOBJ`, 1e-12 curvature), `degenerate_cons`
(duplicate, scaled, dependent, zero and parallel rows; degenerate vertex), `active_set` (many active bounds/rows, strictly and
weakly active multipliers), `bounds_special`, `bounds_huge`, `bounds_tight` (widths 1e-9..1e-3), `scale_q` (1e-6..1e6),
`scale_cols` (column scaling 1e-4..1e4), `scale_rows` (1e-6..1e6), `ill_cond` (cond 1e6..1e12), `maximize`, `obj_const`,
`qform`, `empty_structure` (no rows, empty column, n=1, all fixed, empty rows feasible/infeasible), `portfolio`,
`large_sparse` (n 150..400), `fuzz_mixed`.

For 1460 of the 1780 cases the generator proves the status or the optimal objective itself (82 %), which is what makes `pass` mean more than "agrees with HiGHS".

## Coverage gaps

* No quadratic constraints (the engine rejects them) and no mixed-integer QP (`unsupported` by design): not tested.
* Sizes stop at n = 400; no timing claims (all cases take under 0.13 s in the engine).
* Nonconvex *global* optimality is only checked as "the engine must say `nonconvex`", not against a global solver.
* Ranged rows with an empty range cannot be expressed in MPS (`RANGES` uses `|R|` for L/G rows); such infeasibility is
  tested through crossed column bounds, contradictory parallel rows and inconsistent equalities instead.
* A `UP` bound below zero with the default lower bound is the usual reader special case (lower becomes -inf); the writer
  always emits `UP` before `LO` so crossed bounds survive the file round trip.
* Reproducers are generator cases (smallest of each class), not hand-minimised; `qprepro.py` first checks that the file is
  byte-identical to what the generator writes today.
* HiGHS is the only reference solver; its own defects (above) mean some `pass`es rest on the construction, not on agreement.
