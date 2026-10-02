# Adversarial LP/MILP tests

A seeded, deterministic stress suite for the C++ engine (`src/`) plus a differential harness against HiGHS.
Everything lives in `tools/adv/` (generators, harness, checkers) and `tests/` (runner, seed list, known failures,
reproducers). The engine sources were not touched; bugs found are reported as minimal reproducers, not fixed.

Base commit: `442ca16` (the engine at that commit, built with `g++ -O3 -march=native -std=c++17`). Reference:
HiGHS 1.15.1 through `highspy` (scipy is installed but not used; the false-infeasible presolve behaviour seen with
older scipy builds is handled the same way, see below). All timings come from a shared cloud container (4 cores,
tests run 4 at a time) and are indicative only.

## Running it

```
tests/run_adversarial.sh                # full suite, about 16 minutes on 4 cores
QUICK=1 tests/run_adversarial.sh        # 10 cases per category, about a minute
TARAL=/path/to/engine JOBS=8 OUT=dir tests/run_adversarial.sh
```

It builds the engine, runs every generated case against HiGHS, builds the parser table, replays `tests/repro/` and
prints a pass/fail table by category. Exit status 0 means every failure is a known one (listed in
`tests/known_failures.tsv`, each with a reproducer), 1 means a new failure or a reproducer regression. Needs
`g++`, `python3`, `numpy`, `highspy`. Seeds: `tests/seeds.txt` (case *i* of a category uses seed `SEED_BASE + i`,
`SEED_BASE = 26000 + 1000*k`; the generators use only `random.Random(seed)` and numpy's legacy `RandomState`, and the
generated text was verified identical across `PYTHONHASHSEED` values). Generated files and results are written to
`out/adversarial/` (git-ignored); only generators and seeds are committed.

## What is checked per case

1. The case is written as MPS and the engine is run with `--time-limit 60 --sol --json`.
2. The file is re-read by an **independent oracle parser** (`tools/adv/mpsio.py`, no code shared with `src/mps.cpp`).
   Conventions: RANGES as in the usual readers (E with R<0 gives [rhs+R, rhs]; L/G use |R|; N rows ignored), objective
   constant is minus the RHS on the objective row, `|bound| >= 1e20` is infinite, first named RHS/RANGES/BOUNDS set only.
3. HiGHS solves the same model (default tolerances, presolve on, 1 thread, `mip_rel_gap = 0`,
   `small_matrix_value = 1e-12` because the default silently drops 1e-9-sized coefficients) in a forked child with a hard
   wall-clock cap, since HiGHS 1.15.1 can ignore its time limit or crash on some models.
4. **Gate (all must hold):** same status (optimal / infeasible / unbounded); `|obj_engine - obj_ref| <= 1e-6 * max(1, |obj_ref|)`
   where `obj_engine` is recomputed exactly from the returned point; the returned point violates no row or bound by
   more than `1e-6 * (1 + |violated bound|)` and no integer column by more than 1e-6 (all in exact rational arithmetic
   on the original file, using the `--sol` output). The absolute violation is recorded too; the gate is relative
   because 1e-6 absolute is not attainable in double precision on rows with 1e9 coefficients, whichever solver is used.
5. **Reference false results:** on any disagreement HiGHS is re-solved with presolve off before the engine is blamed.
   A case where the reference is shown wrong is scored `pass*` and says why in `results.jsonl`: presolve-off agrees
   with the engine; the engine's point is exactly feasible and better than the reference objective; the engine matches
   the status proven by construction; or the engine matches an exact rational optimum (near-singular cases). When the
   reference itself hangs, crashes or hits its time limit the case is `other` (unresolved), not a pass.
6. Where the construction proves the status (`expect`), a reference that disagrees with it is not trusted: the engine
   is scored against the construction.

Outcome columns: `pass`, `pass*` (above), `wrong` (wrong status, objective, or an infeasible point), `noans`
(time/iteration limit, `numerical_failure`, `unbounded_relaxation`), `other` (reference unresolved, harness error).
Pass rate = (`pass` + `pass*`) / cases; `noans` and `other` count against it.

## Results (full run from the committed head, 3630 generated cases)

Whole suite: 953 s wall on 4 cores (harness 927 s, parser table and reproducers the rest); two earlier complete runs gave the same counts. 2840 LP and 790 MILP cases; 36 engine runs reached the 60 s cap.

```
category           kind   cases   pass   pass*  wrong  noans  other  known pass-rate
------------------------------------------------------------------------------------
degenerate         LP       200    200       0      0      0      0      0   100.0%
cycling            LP       100    100       0      0      0      0      0   100.0%
dup_coef           LP       150    150       0      0      0      0      0   100.0%
dup_empty_zero     LP       200    200       0      0      0      0      0   100.0%
free_fixed         LP       150    150       0      0      0      0      0   100.0%
huge_bounds        LP       150    127       0     19      4      0     23    84.7%
negzero            LP        60     60       0      0      0      0      0   100.0%
rank_deficient     LP       150    150       0      0      0      0      0   100.0%
redundant_eq       LP       100    100       0      0      0      0      0   100.0%
infeasible         LP       150    150       0      0      0      0      0   100.0%
unbounded          LP       150    150       0      0      0      0      0   100.0%
near_singular      LP       150    100      25     22      3      0     25    83.3%
coef_range         LP       200    160       0     14     26      0     40    80.0%
ranges             LP       200    200       0      0      0      0      0   100.0%
objconst_sense     LP       120    120       0      0      0      0      0   100.0%
neg_lower          LP       100    100       0      0      0      0      0   100.0%
bound_types        LP       100    100       0      0      0      0      0   100.0%
fuzz_mixed         LP       250    250       0      0      0      0      0   100.0%
tolerance_edge     LP       100    100       0      0      0      0      0   100.0%
large_sparse       LP        60     60       0      0      0      0      0   100.0%
knapsack           MILP      80     80       0      0      0      0      0   100.0%
general_int        MILP      80     77       0      0      0      3      3    96.2%
setcover           MILP      60     60       0      0      0      0      0   100.0%
parity             MILP      40     40       0      0      0      0      0   100.0%
fixed_charge       MILP      60     60       0      0      0      0      0   100.0%
int_free_neg       MILP      60     59       0      0      0      1      1    98.3%
markers_mixed      MILP      80     80       0      0      0      0      0   100.0%
int_bound_types    MILP      80     78       2      0      0      0      0   100.0%
objconst_max       MILP      50     50       0      0      0      0      0   100.0%
equality_int       MILP      60     56       0      0      0      4      4    93.3%
unbounded_relax    MILP      40     40       0      0      0      0      0   100.0%
big_values         MILP      40     38       0      0      2      0      2    95.0%
indicator          MILP      60     60       0      0      0      0      0   100.0%
------------------------------------------------------------------------------------
TOTAL              LP      2840   2727      25     55     33      0     88    96.9%
TOTAL              MILP     790    778       2      0      2      8     10    98.7%
pass* = engine matched after the reference was shown wrong/inexact (presolve-off re-solve, exact check, or construction)
wrong = wrong status/objective/feasibility; noans = time/iteration limit, numerical_failure, unresolved MILP; other = harness/generator
```

Reading the table:
* 26 of the 33 generated categories pass completely. Every failure of the engine is in four categories, each with a reproducer:
  `coef_range`, `huge_bounds`, `near_singular` (LP) and `big_values` (MILP). The failures in `general_int`, `int_free_neg` and
  `equality_int` are cases where **neither** solver finishes in 60 s (listed as unresolved, not as an engine/reference disagreement).
* `pass*` is not hidden: in `near_singular` 24 cases count as pass* because the engine matches the exact rational optimum while HiGHS
  (default tolerances) does not, and 1 because the engine's exactly-feasible point has a better objective than the reference's; in
  `int_bound_types` 2 cases count as pass* because HiGHS **presolve** returned a different MIP optimum and presolve-off agrees with
  the engine (`int_bound_types-53045`, `int_bound_types-53054`). Counting those 27 as failures would move the LP rate from 96.9% to
  96.0% and the MILP rate from 98.7% to 98.4%.
* No returned `optimal` point anywhere violated a row or bound by more than 1.6e-7 in absolute terms. The `wrong` column is wrong
  *status* (or objective): in this run the engine never labelled an infeasible point optimal.

## Categories

| category | kind | cases | seeds | what it tests |
|---|---|---|---|---|
| degenerate | LP | 200 | 26000-26199 | assignment and transportation LPs with heavy cost ties, many constraints tight at one vertex (redundant combinations of tight rows) and Klee-Minty cubes: primal/dual degeneracy and stalling |
| cycling | LP | 100 | 27000-27099 | Beale's and Kuhn's classic cycling examples (row/column permuted, power-of-two scaled) and homogeneous degenerate LPs that make textbook pivot rules cycle |
| dup_coef | LP | 150 | 28000-28149 | the same (row, column) coefficient listed two or three times in COLUMNS (sums, pairs that cancel to zero, duplicate objective entries) |
| dup_empty_zero | LP | 200 | 29000-29199 | duplicate / scaled / sign-flipped rows, empty rows (satisfied and violated), explicit zero entries, zero rows, empty and zero columns, duplicate columns |
| free_fixed | LP | 150 | 30000-30149 | free (FR) and fixed (FX or LO=UP) variables, incl. free columns with no entries and fixed columns that shift row activities |
| huge_bounds | LP | 150 | 31000-31149 | bounds, RHS and RANGES written as 1e20/1e30 (infinity by common convention), huge finite bounds (1e15), an unbounded ray bounded only by the sentinel |
| negzero | LP | 60 | 32000-32059 | '-0' spelled in coefficients, RHS, RANGES, bounds, the objective constant and costs, incl. UP -0 on a column with the default lower bound |
| rank_deficient | LP | 150 | 33000-33149 | equality systems whose rank is below the row count: consistent dependent rows (feasible) |
| redundant_eq | LP | 100 | 34000-34099 | redundant equalities (consistent, feasible) and inconsistent copies (infeasible) |
| infeasible | LP | 150 | 35000-35149 | structured infeasible LPs: contradictory bounds, Farkas rows, equalities without nonnegative solution, a row beyond the bound box, infeasible plus an improving ray, disjoint ranges on identical rows |
| unbounded | LP | 150 | 36000-36149 | unbounded LPs with a recession ray (row-compatible signs, free column), MIN and MAX |
| near_singular | LP | 150 | 37000-37149 | square equality systems with condition number 1e8-1e12 (prescribed singular values, Hilbert, nearly dependent rows), plain and objective-stable costs; an exact rational optimum arbitrates |
| coef_range | LP | 200 | 38000-38199 | exactly rescaled models (row scales 10^-4..10^4, column scales 10^-5..10^5) so entries span 1e-9..1e9 while the optimum is unchanged |
| ranges | LP | 200 | 39000-39199 | RANGES on E/L/G rows (positive, negative and zero R) with mostly negative RHS, plus an N row carrying RHS and RANGES |
| objconst_sense | LP | 120 | 40000-40119 | objective constants (RHS on the objective row: +-7, 1e6, 0.001, ...) and OBJSENSE MAX/MIN in section, inline and MAXIMIZE spellings |
| neg_lower | LP | 100 | 41000-41099 | negative lower bounds, negative upper bounds (MI+UP), boxes straddling zero, free and fixed columns |
| bound_types | LP | 100 | 42000-42099 | LO/UP/FX/FR/MI/PL combinations and orders on LP columns |
| fuzz_mixed | LP | 250 | 43000-43249 | random mixtures of the above (up to ~30x40), duplicates, zeros, ranges, constants, MAX |
| tolerance_edge | LP | 100 | 44000-44099 | infeasible by exactly 1e-2..1e-5 (a row demanding the box maximum plus delta) and the delta=0 feasible twin |
| large_sparse | LP | 60 | 45000-45059 | sparse LPs with 60-150 rows and 80-200 columns mixing every row and column kind |
| knapsack | MILP | 80 | 46000-46079 | binary knapsacks (BV bounds) with 1-3 rows, up to 40 items |
| general_int | MILP | 80 | 47000-47079 | general integers with box/lower bounds via MARKER blocks or LI/UI bound types |
| setcover | MILP | 60 | 48000-48059 | set covering with binary columns, up to 45 sets x 25 elements |
| parity | MILP | 40 | 49000-49039 | even-coefficient equalities with odd RHS: LP-feasible but integer-infeasible (and feasible twins) |
| fixed_charge | MILP | 60 | 50000-50059 | fixed-charge networks with big-M links x <= M y, M up to 1e5 |
| int_free_neg | MILP | 60 | 51000-51059 | integers with free / negative / upper-negative bounds made finite |
| markers_mixed | MILP | 80 | 52000-52079 | mixed integer/continuous columns, interleaved MARKER blocks, duplicate coefficients, BV/LI/UI versus markers |
| int_bound_types | MILP | 80 | 53000-53079 | BV, LI, UI, MI+UI and fractional LI/UI values (ceil/floor semantics) on integer columns |
| objconst_max | MILP | 50 | 54000-54049 | MILP with objective constants and OBJSENSE MAX |
| equality_int | MILP | 60 | 55000-55059 | integer equality systems |
| unbounded_relax | MILP | 40 | 56000-56039 | unbounded LP relaxation: integer-feasible ray (unbounded) versus no integer point (infeasible) |
| big_values | MILP | 40 | 57000-57039 | MILP with all coefficients, bounds and costs scaled by 1e3..1e7 |
| indicator | MILP | 60 | 58000-58059 | big-M indicator links with M up to 1e7 and binary y: integrality tolerance against big-M |

Plus the parser strictness table (88 inputs) and 19 reproducers. Generated: 2840 LP + 790 MILP (the brief asked for at least 2000 and 500).

## MPS parser strictness: engine vs HiGHS reader

`tools/adv/parser_table.py` feeds each input to the engine and to HiGHS' reader (then solves what HiGHS read, presolve off).
"Spec says" is **our reading of the MPS format** (accept / reject / either). Where the format is silent or readers legitimately
differ the verdict is `either`, the input counts as pass, and the table shows the difference. "Who is right" follows from that
reading plus the expected optimum where one exists. Several rows are judgement calls (the `why` strings in the script say which);
treat the table as evidence, not a ruling. Result: 88 inputs, 80 pass, 8 fail.

The known case, `BOUNDS` on an undeclared column: **the engine rejects (parse_error), HiGHS ignores the line and solves.** The format
does not spell this case out, but a bound refers to a column and none was declared, and other readers report an error, so we score
the engine as right (strict) and HiGHS as lenient, with that caveat. The eight failures are the engine being too lenient or wrong
in its output: `missing_endata`, `dup_row_name`, `dup_column_noncontiguous` (the engine merges the two blocks, HiGHS makes a second
column, so they solve different models), `dup_section_columns`, `num_hex_value`, and three non-finite coefficient inputs
(`nan`, `inf`, `1e400`) where the engine prints `optimal` with a NaN objective and writes **invalid JSON** (a bare `nan`).

| input | what | spec says | ours | HiGHS | who is right | ours verdict |
|---|---|---|---|---|---|---|
| `free_control` | control: plain free format | accept | optimal 3 | optimal 3 | both | PASS |
| `fixed_aligned` | strict fixed columns, names <= 8 chars | accept | optimal 3 | optimal 3 | both | PASS |
| `fixed_names_with_spaces` | fixed columns, names containing spaces | accept | optimal 3 | optimal 3 | both | PASS |
| `fixed_shifted_fields` | fixed-format file whose fields are shifted one column but whitespace separated | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `rhs_without_set_name` | RHS line with the set name omitted (row value) | accept | optimal 3 | optimal 3 | both | PASS |
| `rhs_two_pairs_one_line` | RHS with two row/value pairs on one line | accept | optimal 3 | optimal 3 | both | PASS |
| `lower_case_row_type` | row type in lower case ('g') | either | optimal 3 | error | spec silent: both defensible | PASS |
| `missing_endata` | no ENDATA line | reject | optimal 3 | error | HiGHS only | FAIL |
| `text_after_endata` | garbage after ENDATA | accept | optimal 3 | optimal 3 | both | PASS |
| `missing_name_line` | no NAME line | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `empty_file` | empty file | reject | error | error | both | PASS |
| `whitespace_only_file` | only blank lines | reject | error | error | both | PASS |
| `no_columns_section` | no COLUMNS section | either | optimal 0 | optimal 0 | same behaviour | PASS |
| `no_rhs_section` | no RHS section (all zero) | accept | optimal 0 | optimal 0 | both | PASS |
| `bounds_before_columns` | BOUNDS section before COLUMNS | reject | error | optimal 3 | ours only | PASS |
| `unknown_section` | unknown section name | reject | error | optimal 0 | ours only | PASS |
| `no_objective_row` | no N row at all | either | error | optimal 0 | spec silent: both defensible | PASS |
| `second_objective_row` | two N rows: the first is the objective, the second is a free row | accept | optimal 3 | optimal 3 | both | PASS |
| `rows_declared_unused` | declared row with no entries | accept | optimal 3 | optimal 3 | both | PASS |
| `dup_row_name` | two rows with the same name | reject | optimal 3 | optimal 3 | neither | FAIL |
| `dup_column_noncontiguous` | column x listed, then y, then x again | reject | optimal 1.5 | optimal 0 | neither | FAIL |
| `dup_entry_same_cell` | the same (row, column) listed twice | either | optimal 1.5 | optimal 3 | spec silent: both defensible | PASS |
| `dup_rhs_entry` | RHS for the same row twice | either | optimal 5 | optimal 3 | spec silent: both defensible | PASS |
| `dup_bound_same_type` | UP bound listed twice on one column | either | optimal 4 | optimal 3 | spec silent: both defensible | PASS |
| `dup_section_columns` | COLUMNS section twice | reject | optimal 3 | optimal 3 | neither | FAIL |
| `row_and_column_same_name` | a row and a column both named x | accept | optimal 3 | optimal 3 | both | PASS |
| `name_equals_section_keyword` | a column named RHS and a row named BOUNDS | accept | optimal 3 | optimal 3 | both | PASS |
| `tabs_free` | tab-separated fields | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `tab_leading` | data lines start with a tab | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `crlf` | CRLF line endings | accept | optimal 3 | optimal 3 | both | PASS |
| `cr_only` | CR-only line endings | either | error | error | same behaviour | PASS |
| `no_final_newline` | no newline after ENDATA | accept | optimal 3 | optimal 3 | both | PASS |
| `trailing_spaces` | trailing blanks on every line | accept | optimal 3 | optimal 3 | both | PASS |
| `utf8_bom` | UTF-8 byte-order mark before NAME | either | error | optimal 3 | spec silent: both defensible | PASS |
| `blank_lines_inside` | blank lines between and inside sections | accept | optimal 3 | optimal 3 | both | PASS |
| `full_line_comments` | '*' comment lines everywhere | accept | optimal 3 | optimal 3 | both | PASS |
| `indented_star_comment` | comment with leading space ' * text' inside COLUMNS | reject | error | error | both | PASS |
| `dollar_inline_comment` | '$' inline comment after the data | either | error | optimal 3 | spec silent: both defensible | PASS |
| `long_names_40` | 40-character names in free format | accept | optimal 3 | optimal 3 | both | PASS |
| `long_names_300` | 300-character names in free format | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `names_special_chars` | names such as x[1], a.b-c, 7up | accept | optimal 3 | optimal 3 | both | PASS |
| `very_long_line` | one data line with 100000 spaces of padding | accept | optimal 3 | optimal 3 | both | PASS |
| `num_fortran_D_exponent` | coefficient written as '1D0' | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `num_plus_sign` | coefficient written as '+1' | accept | optimal 3 | optimal 3 | both | PASS |
| `num_leading_dot` | coefficient written as '.5E1' | accept | optimal 6 | optimal 6 | both | PASS |
| `num_trailing_dot` | coefficient written as '1.' | accept | optimal 3 | optimal 3 | both | PASS |
| `num_empty_exponent` | coefficient written as '1e' | reject | error | optimal 3 | ours only | PASS |
| `num_nan_value` | coefficient written as 'nan' | reject | optimal nan | optimal nan | n/a (engine reports a non-finite or unreadable result whatever the input convention) | FAIL |
| `num_hex_value` | coefficient written as '0x10' | reject | optimal 6 | optimal 6 | neither | FAIL |
| `num_decimal_comma` | coefficient written as '1,5' | reject | error | optimal 3 | ours only | PASS |
| `num_trailing_garbage` | coefficient written as '3abc' | reject | error | optimal 6 | ours only | PASS |
| `num_double_sign` | coefficient written as '--3' | reject | error | optimal 0 | ours only | PASS |
| `num_inf_value` | coefficient written as 'inf' | either | optimal nan | optimal 6 | n/a (engine reports a non-finite or unreadable result whatever the input convention) | FAIL |
| `num_overflow_1e400` | coefficient written as '1e400' | either | optimal nan | optimal 6 | n/a (engine reports a non-finite or unreadable result whatever the input convention) | FAIL |
| `bound_undeclared_column` | BOUNDS line for a column that is not in COLUMNS (known case) | reject | error | optimal 3 | ours only | PASS |
| `bound_unknown_type` | bound type XX | reject | error | error | both | PASS |
| `bound_lower_case_type` | bound type 'up' | either | error | error | same behaviour | PASS |
| `bound_up_negative_default_lo` | UP -1 on a column with default lower bound 0 | either | optimal 7 | infeasible | spec silent: both defensible | PASS |
| `bound_mi_default_upper` | MI on a column (upper bound afterwards?) | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `bound_bv_with_value` | BV with an (ignored) value | accept | optimal 5 | optimal 5 | both | PASS |
| `bound_lo_gt_up` | LO 5 and UP 3 on one column | accept | infeasible | infeasible | both | PASS |
| `bound_two_sets` | two named BOUNDS sets (second sets x <= 1) | either | optimal 3 | optimal 3 | same behaviour | PASS |
| `bound_sc_semicontinuous` | SC (semi-continuous) bound | either | error | optimal 3 | spec silent: both defensible | PASS |
| `bound_free_then_up` | FR followed by UP 2 on x | accept | optimal 4 | optimal 3 | ours only | PASS |
| `integer_marker_no_bounds` | integer column in a MARKER block with no bound, minimising -x | either | unbounded | optimal -1 | spec silent: both defensible | PASS |
| `marker_unbalanced` | INTORG marker never closed | either | optimal 3 | infeasible | spec silent: both defensible | PASS |
| `marker_bad_token` | MARKER line with an unknown type | reject | error | error | both | PASS |
| `bound_missing_value` | UP bound with no value | reject | error | error | both | PASS |
| `bound_bad_number` | UP bound with a non-number | reject | error | optimal 6 | ours only | PASS |
| `rhs_missing_value` | RHS entry with a row but no value | reject | error | error | both | PASS |
| `ranges_unknown_row` | RANGES for a row that does not exist | reject | error | optimal 3 | ours only | PASS |
| `rows_line_three_fields` | ROWS line with a third field | reject | error | optimal 0 | ours only | PASS |
| `empty_optional_sections` | empty RHS, RANGES and BOUNDS sections | accept | optimal 0 | optimal 0 | both | PASS |
| `name_line_extra_words` | NAME line with several words | accept | optimal 3 | optimal 3 | both | PASS |
| `rhs_unknown_row` | RHS for a row that does not exist | reject | error | optimal 3 | ours only | PASS |
| `columns_unknown_row` | COLUMNS entry for an undeclared row | reject | error | optimal 3 | ours only | PASS |
| `range_negative_on_E` | E row with R<0: interval [rhs+R, rhs] | accept | optimal 3 | optimal 3 | both | PASS |
| `range_on_G_negative` | G row with R<0: /R/ is used, interval [rhs, rhs+/R/] | accept | optimal -7 | optimal -7 | both | PASS |
| `range_on_N_row` | RANGES entry for the objective row (ignored) | accept | optimal 3 | optimal 3 | both | PASS |
| `range_zero` | R = 0 on an L row (becomes equality) | accept | optimal -5 | optimal -5 | both | PASS |
| `rhs_on_objective_row` | RHS on the objective row (constant term, sign: obj = c'x - rhs) | accept | optimal -7 | optimal -7 | both | PASS |
| `objsense_max_section` | OBJSENSE section with MAX on its own line | accept | optimal 4 | optimal 4 | both | PASS |
| `objsense_inline` | 'OBJSENSE MAX' on one line | either | optimal 4 | optimal 4 | same behaviour | PASS |
| `objsense_maximize` | OBJSENSE with the word MAXIMIZE | accept | optimal 4 | optimal 4 | both | PASS |
| `objsense_lowercase` | OBJSENSE with lower-case 'max' | either | error | optimal 4 | spec silent: both defensible | PASS |
| `objsense_garbage` | OBJSENSE with an unknown word | reject | error | optimal 0 | ours only | PASS |
| `columns_line_four_tokens` | COLUMNS line with an odd number of fields | reject | error | optimal 6 | ours only | PASS |
| `columns_line_seven_tokens` | COLUMNS line with three pairs | reject | error | optimal 3 | ours only | PASS |

## Known failures

Each is listed per case in `tests/known_failures.tsv` and isolated in `tests/repro/` (replayed by the run script; expected status
comes from the reference or from an exact derivation, never from the engine). Notes on the other engine methods come from running
`--method dual` and `--method ipm` on the reproducers. Nothing in `src/` was changed and no root cause was established; the
paragraphs describe symptoms and likely directions only.

### coef_range_stall (25 cases)

Reproducer `tests/repro/coef_range_stall.mps`: expected optimal -32; observed: default simplex runs to the time limit (tens of millions of iterations in a few seconds) on a 3-row, 4-column LP. Cases: coef_range-38013, coef_range-38019, coef_range-38078, coef_range-38084, coef_range-38090, coef_range-38098, ....

A 3 x 4 LP whose coefficients span 1e-4..3e4 (an exactly rescaled version of a well-behaved model). The default primal simplex spends the whole time limit (about 8e6 iterations per second, i.e. it is cycling or stalling, not slowly converging) while HiGHS and the engine's own --method dual (4 iterations, objective -32) and --method ipm solve it instantly. This is the dominant failure mode in the coef_range category (see docs/adversarial.md for counts). Not root-caused.

### coef_range_false_infeasible (11 cases)

Reproducer `tests/repro/coef_range_false_infeasible.mps`: expected optimal -70; observed: default simplex returns infeasible after 2 iterations. Cases: coef_range-38007, coef_range-38039, coef_range-38057, coef_range-38069, coef_range-38079, coef_range-38095, ....

One column and two rows with coefficients 20 and 4e9, RHS -1e-3 and -3.2e5, and upper bound 7e-5. The rows give x >= 5e-5 and x <= 8e-5, the bound gives x <= 7e-5, so the feasible set is [5e-5, 7e-5] and the optimum of min -1e6 x is x = 7e-5, objective -70. HiGHS (simplex and interior point, presolve on and off) returns -70, and so do the engine's --method dual (-70) and --method ipm (-69.9999999916), but the default primal simplex declares the model infeasible after 2 iterations. Not root-caused: the variable lives at the 1e-5 scale while its row coefficient is 4e9, which points at an absolute feasibility tolerance applied to an unscaled model.

### coef_range_false_unbounded (3 cases)

Reproducer `tests/repro/coef_range_false_unbounded.mps`: expected optimal 0; observed: default simplex returns unbounded after 1 iteration. Cases: coef_range-38040, coef_range-38049, coef_range-38071.

One row 0.1 C0 + 2e8 C1 >= 2e5 with C0, C1 >= 0 and objective min 1e4 C1. The cost is non-negative and C1 = 0, C0 = 2e6 is feasible, so the optimum is 0 and the LP is bounded. The default primal simplex returns 'unbounded' after one iteration; --method dual returns optimal 0 and so does HiGHS. A wrong unbounded claim is worse than a failure to finish because it is a definite (false) statement about the model. Not root-caused: the row coefficients differ by nine orders of magnitude.

### coef_range_numfail (1 cases)

Reproducer `tests/repro/coef_range_numfail.mps`: expected optimal 0; observed: default simplex returns numerical_failure: phase 1 direction without a blocking variable. Cases: coef_range-38127.

A 5 x 5 LP with equality rows and coefficients 5e-5..5e9 that has optimal objective 0. The default primal simplex stops with numerical_failure ('phase 1 direction without a blocking variable') after 5 iterations; --method dual and --method ipm return optimal 0, as does HiGHS. Not root-caused.

### huge_bounds_sentinel_unbounded (19 cases)

Reproducer `tests/repro/huge_bounds_sentinel_unbounded.mps`: expected unbounded; observed: optimal, objective -5e+30. Cases: huge_bounds-31000, huge_bounds-31008, huge_bounds-31016, huge_bounds-31024, huge_bounds-31032, huge_bounds-31040, ....

No rows, one column with cost 5 and bounds [-1e30, 1e30]. By the common convention (the common readers, HiGHS included, all treat |bound| >= 1e20 as infinity) the column is free, so the LP is unbounded below. The engine's simplex and dual paths treat 1e30 as a finite bound and report 'optimal' with objective -5e+30, a point nobody can use; only --method ipm reports unbounded. The MPS format itself does not define the sentinel, so this is a convention gap rather than a format violation, but returning 'optimal' with a 1e30 objective is a wrong answer under every common reading.

### huge_bounds_sentinel_numfail (4 cases)

Reproducer `tests/repro/huge_bounds_sentinel_numfail.mps`: expected optimal -4; observed: numerical_failure: final point violates constraints by 0.8889. Cases: huge_bounds-31019, huge_bounds-31068, huge_bounds-31089, huge_bounds-31095.

A 4 x 4 LP in which three columns are bounded by +-1e20 and two rows carry RANGES 1e20 (so one side of each is infinite by the 1e20 convention). HiGHS and the engine's --method dual / --method ipm return optimal -4. The default simplex treats the 1e20 values as finite data, runs 5 iterations and stops with numerical_failure ('final point violates constraints by 0.8889'): the huge finite bounds enter the arithmetic of the ratio test and destroy the basic solution. Same convention gap as the other huge_bounds reproducers; with 1e30 in place of 1e20 the failure is identical.

### near_singular_false_infeasible (6 cases)

Reproducer `tests/repro/near_singular_false_infeasible.mps`: expected optimal 10; observed: default simplex and --method dual return infeasible; --method ipm returns optimal 9.99999998998. Cases: near_singular-37005, near_singular-37083, near_singular-37097, near_singular-37103, near_singular-37110, near_singular-37119.

A 4 x 4 square equality system B x = b, x >= 0, with B of condition number about 1e11 (prescribed singular values) written with full 17-digit coefficients. In exact rational arithmetic the solution is unique and nonnegative, so the LP is feasible with objective 9.999999989977592. HiGHS (default tolerances) returns optimal 9.99999998251 (a point within its 1e-9 feasibility tolerance), the engine's --method ipm returns 9.99999998998, but the default simplex and the dual method declare the model infeasible. At this conditioning the feasibility verdict is tolerance-dependent; the reproducer records that the engine answers 'infeasible' for a model that is exactly feasible and that other solvers accept.

### near_singular_objective (16 cases)

Reproducer `tests/repro/near_singular_objective.mps`: expected optimal 38; observed: optimal, objective 36.2183747274 instead of 38.0. Cases: near_singular-37003, near_singular-37027, near_singular-37045, near_singular-37056, near_singular-37057, near_singular-37062, ....

A 4 x 4 square equality system with condition number about 1e11, so the feasible set is a single point and the optimal objective is exactly 38. The engine returns optimal 36.2183747274 (relative error 0.047) and HiGHS returns 36.2172385 (also off by 0.047): both stop inside the cloud of points that satisfy the equalities to their 1e-9 feasibility tolerance. The engine's answer is infeasible by about 1e-9 relative to the data but is nowhere near the exact optimum, and it differs from HiGHS by 3e-5, above the 1e-6 differential gate. This reproducer documents an inherent limit of tolerance-based LP solvers at this conditioning, not a unique defect of the engine.

### near_singular_stall (3 cases)

Reproducer `tests/repro/near_singular_stall.mps`: expected optimal 15; observed: default simplex runs to the time limit (millions of iterations) on a 4x4 LP. Cases: near_singular-37068, near_singular-37098, near_singular-37139.

A 4 x 4 square equality system with condition number about 1e10 and a unique feasible point (objective 14.999999242). HiGHS solves it; the engine's default simplex iterates until the time limit (over 2 million iterations in 3 seconds), i.e. it cycles or stalls on a model that needs four pivots.

### big_values_numfail (2 cases)

Reproducer `tests/repro/big_values_numfail.mps`: expected optimal -8.782e+09; observed: numerical_failure after 14 nodes although its incumbent (-8.782e9) is already the exact optimum. Cases: big_values-57025, big_values-57027.

A 5-integer, 1-continuous MIP with coefficients and costs of order 1e7 (objective about -8.8e9). Exhaustive enumeration in exact rational arithmetic gives the optimum -8.782e9 (C = 6,5,24,1,16 and C5 = 10.4) and HiGHS at default options returns the same. The engine finds exactly that incumbent but ends with numerical_failure ('node LP numerical_failure: final point violates constraints by 0.000000') and a best bound of -8.851e9, so it cannot prove optimality; the same happens with --method dual and --method ipm because the node LP path is shared. Side finding about the reference: HiGHS 1.15.1 with mip_feasibility_tolerance set to 1e-9 returns -8.256e9 (and -8.678e9 with some other options), both worse than the true optimum, which is why the harness leaves HiGHS tolerances at their defaults. Engine side not root-caused.

### Same convention gap, 1e20 spelling
`tests/repro/huge_bounds_sentinel_1e20_unbounded.mps` (hand-written): see `huge_bounds_sentinel_unbounded`.

### Parser reproducers
`tests/repro/parser_*.mps` (8 files), expected status `parse_error`; see the parser section and the `report` field of each entry in
`tests/repro/expected.json`.

### Unresolved: neither solver finishes (not an engine/reference disagreement)
`general_int` (3), `int_free_neg` (1), `equality_int` (4): both the engine and HiGHS hit the 60 s limit on integer equality and
general-integer models up to about 30 x 40. They are class `both_time_limit` in `known_failures.tsv` and count against the pass rate.
No reproducer is possible because there is no reference answer.

## Reference-solver observations (HiGHS 1.15.1)
* Presolve returned a different MIP optimum on 2 `int_bound_types` cases (presolve off agrees with the engine).
* With `mip_feasibility_tolerance = 1e-9` HiGHS returned -8.256e9 (and -8.678e9 with other options) on the `big_values_numfail` model,
  whose exact optimum is -8.782e9; default tolerances are correct. In that tighter configuration HiGHS also ignored its time limit
  (hang) on `big_values-57014` and reported `kSolveError` on `big_values-57027`. The harness therefore keeps default tolerances and runs
  the reference in a forked child with a hard wall cap.
* The default `small_matrix_value = 1e-9` silently drops 1e-9-sized coefficients, which turned exact `coef_range` models into different
  models; the harness sets it to 1e-12.

## What is not covered
* **QP and the interior-point path as a primary method are not tested, and `--method dual` / `--method ipm` are not run through the
  differential harness** (they were only run on the reproducers). The harness drives the default CLI path: primal simplex for LP,
  branch and bound for MILP. Many `coef_range` and `huge_bounds` reproducers are solved by the dual and interior-point paths, so a
  rerun of the same seeds with those methods would change the picture; it was not done.
* Sizes are small: LP up to 150 rows x 200 columns, MILP mostly up to about 45 columns. No industrial-scale models, no performance
  or memory claims. Timings come from a shared cloud container, four cases at a time, and are indicative only.
* Generated cases are free-format; fixed-column format appears only in a handful of parser-table inputs.
* Not generated: SOS sets, semicontinuous variables (one parser input only), indicator or quadratic sections, `--node-limit`,
  warm starts, time-limit behaviour beyond observing the 60 s cap, files above a few hundred KB, non-text input beyond one BOM case.
* The checker validates **primal feasibility and the objective only**. An `optimal` claim is not certified by a dual solution; an
  `infeasible` or `unbounded` claim is checked only by agreement with the reference (re-solved without presolve) or by construction.
  Where HiGHS is itself inexact (`near_singular`) the comparison is arbitrated by an exact rational optimum, and those cases are
  tolerance-defined: above condition number about 1e10 any tolerance-based solver, HiGHS included, ends in the cloud of nearly feasible
  points. The 25 `near_singular` failures therefore mean "differs from the exact optimum or from HiGHS by more than 1e-6", not that the
  engine is uniquely wrong.
* The violation gate is relative to `1 + |bound|`; an absolute 1e-6 gate is not attainable in double precision on rows with 1e9
  coefficients for any solver.
* HiGHS 1.15.1 via `highspy` is the only reference. scipy's bundled HiGHS (the build with the presolve false-infeasible problem mentioned
  in the brief) is not used. When a model is both infeasible and unbounded the status is settled by presolve-off reference runs only.
* Determinism: case generation is reproducible bit for bit. Engine results depend on compiler flags (`-march=native`) and can move
  borderline `near_singular` cases between machines.
* One run per case: no repeats, no sanitizer builds, no fuzzing of malformed binary input.
