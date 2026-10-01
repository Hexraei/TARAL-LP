#!/usr/bin/env python3
"""MPS semantics regression suite for the C++ engine.

Each case is a tiny hand-written MPS model that isolates one convention (free-variable names that
share a 7-character prefix, ranges on every row type, negative upper bounds, the objective-row RHS
constant, multiple RHS/BOUNDS sets, OBJSENSE, fixed-column names containing spaces, integer markers,
quadratic objectives). Every case states the answer derived by hand; when highspy is installed the
same original file is also solved by HiGHS and must agree. HiGHS is a reference only.

Usage: python benchmarks/mps_semantics_check.py [--engine ./taral]   (exit code 1 on any failure)
"""
import argparse, json, os, subprocess, sys, tempfile


def fixed(kind="", f1="", f2="", f3="", f4="", f5=""):
    """One fixed-column MPS data line: fields at columns 2-3, 5-12, 15-22, 25-36, 40-47, 50-61."""
    line = " " + kind.ljust(2) + " " + f1.ljust(8) + "  " + f2.ljust(8) + "  " + f3.rjust(12)
    if f4:
        line += "   " + f4.ljust(8) + "  " + f5.rjust(12)
    return line.rstrip()


# Documented convention differences: HiGHS (1.15.1) reads these files differently, so its answer is reported
# but not required to match. The engine follows the convention used by the independent original-model checker.
CONVENTION_DIFFERS = {
    "negative_upper_bound": "HiGHS keeps the default lower bound 0 (infeasible); engine and checker use -inf",
}

# name -> (mps text, expected status, expected objective or None, why)
CASES = {
    "free_names_shared_prefix": ("""NAME FRNAMES
ROWS
 N obj
 E r1
 E r2
 L r3
COLUMNS
 ABCDEFG1 obj 1 r1 1
 ABCDEFG2 obj 2 r2 1
 ABCDEFGP obj 1 r3 1
RHS
 rhs r1 1 r2 2
 rhs r3 3
BOUNDS
 FR bnd ABCDEFG1
 FR bnd ABCDEFG2
ENDATA
""", "optimal", 5.0, "two free columns sharing a 7-char prefix plus a column named like a split clone; x1=1, x2=2, ABCDEFGP=0"),
    "ranges_all_row_types": ("""NAME RANGES
ROWS
 N obj
 E e_pos
 E e_neg
 L l_rng
 G g_rng
COLUMNS
 a obj 1 e_pos 1
 b obj -1 e_neg 1
 c obj 1 l_rng 1
 d obj -1 g_rng 1
RHS
 rhs e_pos 2 e_neg 5
 rhs l_rng 10 g_rng 1
RANGES
 rng e_pos 3 e_neg -4
 rng l_rng 4 g_rng 6
ENDATA
""", "optimal", 2 - 5 + 6 - 7, "E+R: [2,5] min a=2; E-R: [1,5] max b=5; L: [6,10] min c=6; G: [1,7] max d=7"),
    "negative_upper_bound": ("""NAME NEGUP
ROWS
 N obj
 L r
COLUMNS
 x obj -1 r 1
RHS
 rhs r 100
BOUNDS
 UP bnd x -2
ENDATA
""", "optimal", 2.0, "UP < 0 with default lower 0 makes the lower bound -inf; max x = -2"),
    "objective_rhs_constant": ("""NAME OBJCONST
ROWS
 N obj
 G r
COLUMNS
 x obj 1 r 1
RHS
 rhs r 3 obj -7
ENDATA
""", "optimal", 3 + 7.0, "RHS -7 on the objective row adds +7"),
    "first_rhs_and_bound_set_only": ("""NAME SETS
ROWS
 N obj
 G r
COLUMNS
 x obj 1 r 1
RHS
 rhs1 r 4
 rhs2 r 9
BOUNDS
 LO bnd1 x 1
 LO bnd2 x 6
ENDATA
""", "optimal", 4.0, "only the first RHS set (4) and first BOUNDS set (lo 1) apply"),
    "objsense_max_next_line": ("""NAME MAXNL
OBJSENSE
    MAX
ROWS
 N obj
 L r
COLUMNS
 x obj 3 r 1
RHS
 rhs r 2
ENDATA
""", "optimal", 6.0, "maximise 3x with x <= 2"),
    "objsense_max_same_line": ("""NAME MAXSL
OBJSENSE MAX
ROWS
 N obj
 L r
COLUMNS
 x obj 3 r 1
RHS
 rhs r 2 obj -1
ENDATA
""", "optimal", 7.0, "maximise 3x + 1 with x <= 2"),
    "fixed_columns_names_with_spaces": ("\n".join([
        "NAME          SPACES",
        "ROWS",
        fixed("N", "COST"),
        fixed("G", "ROW A"),
        fixed("L", "ROW B"),
        "COLUMNS",
        fixed("", "X ONE", "COST", "1", "ROW A", "1"),
        fixed("", "X TWO", "COST", "2", "ROW A", "1"),
        fixed("", "X TWO", "ROW B", "1"),
        "RHS",
        fixed("", "RHS", "ROW A", "3", "ROW B", "5"),
        "BOUNDS",
        fixed("UP", "BND", "X ONE", "1"),
        "ENDATA", ""]), "optimal", 1 + 2 * 2.0, "names with spaces need fixed columns; x1=1 (bounded), x2=2"),
    "bound_types_mi_fx_pl": ("""NAME BOUNDS
ROWS
 N obj
 G r1
 G r2
COLUMNS
 a obj 1 r1 1
 b obj 1 r2 1
 c obj 1
RHS
 rhs r1 -4 r2 -1
BOUNDS
 MI bnd a
 FX bnd c 2.5
 PL bnd b
ENDATA
""", "optimal", -4 + 0 + 2.5, "MI frees the lower bound (a=-4), b>=0 so b=0, FX c=2.5"),
    "free_format_d_exponent": ("""NAME DEXP
ROWS
 N obj
 G r
COLUMNS
 x obj 1.5D0 r 2.0D+00
RHS
 rhs r 1.0D1
ENDATA
""", "optimal", 7.5, "D exponents; 2x >= 10 so x=5, cost 1.5*5"),
    "integer_marker_rejected_by_lp": ("""NAME INTS
ROWS
 N obj
 G r
COLUMNS
 M1 'MARKER' 'INTORG'
 x obj 1 r 2
 M2 'MARKER' 'INTEND'
RHS
 rhs r 3
BOUNDS
 UP bnd x 10
ENDATA
""", "unsupported", None, "integer columns are never silently relaxed by the LP engine (MIP optimum would be 2)"),
    "quadratic_objective_rejected_by_lp": ("""NAME QUAD
ROWS
 N obj
 G r
COLUMNS
 x obj 1 r 1
RHS
 rhs r 1
QUADOBJ
 x x 2
ENDATA
""", "unsupported", None, "a quadratic objective is never silently dropped by the LP engine"),
}


def highs_objective(path):
    try:
        import highspy
    except ImportError:
        return None
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(path)
    h.run()
    status = h.modelStatusToString(h.getModelStatus())
    return status, h.getInfo().objective_function_value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="./taral")
    a = ap.parse_args()
    failures = 0
    with tempfile.TemporaryDirectory() as d:
        for name, (text, want_status, want_obj, why) in CASES.items():
            mps, js = os.path.join(d, name + ".mps"), os.path.join(d, name + ".json")
            open(mps, "w").write(text)
            subprocess.run([a.engine, mps, "--time-limit", "10", "--json", js], capture_output=True)
            got = json.load(open(js))
            ok = got["status"] == want_status
            if want_obj is not None:
                ok = ok and got["objective"] is not None and abs(got["objective"] - want_obj) <= 1e-9 * (1 + abs(want_obj))
            ref = highs_objective(mps)
            ref_note = ""
            if ref is not None:
                hs, ho = ref
                if name in CONVENTION_DIFFERS:
                    ref_note = f" | HiGHS {hs} {ho:.10g} (convention difference: {CONVENTION_DIFFERS[name]})"
                elif want_obj is not None:
                    ref_ok = hs == "Optimal" and abs(ho - want_obj) <= 1e-9 * (1 + abs(want_obj))
                    ok = ok and ref_ok
                    ref_note = f" | HiGHS {hs} {ho:.10g}" + ("" if ref_ok else " DISAGREES")
                else:
                    ref_note = f" | HiGHS {hs} {ho:.10g} (reference only)"
            failures += not ok
            print(f"{'PASS' if ok else 'FAIL'} {name}: engine {got['status']} {got['objective']} "
                  f"(expected {want_status} {want_obj}){ref_note}  -- {why}")
    print(f"{len(CASES) - failures}/{len(CASES)} semantic cases pass")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
