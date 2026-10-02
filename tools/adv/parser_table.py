#!/usr/bin/env python3
"""MPS parser strictness table: the same input fed to the engine and to HiGHS, with a spec verdict per input.

Each case states what the MPS format (IBM/CPLEX description as implemented by the common readers) requires:
  accept  -> a conforming reader must read it and produce `expect` (status, objective)
  reject  -> a conforming reader must report an error (the input is malformed)
  either  -> the format is silent or readers legitimately differ; any behaviour is defensible, differences are listed
The verdict column is *our reading of the format*; where it is a judgement call the `why` text says so.
Usage: parser_table.py --engine ./taral --out DIR   (writes DIR/parser_table.md and DIR/parser_table.json)
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref import read_ref, solve_ref
from mpsio import parse_mps, OracleError

BASE = """NAME T
ROWS
 N obj
 G r1
COLUMNS
 x obj 1 r1 1
 y obj 2 r1 1
RHS
 rhs r1 3
ENDATA
"""
OPT3 = ("optimal", 3.0)  # min x+2y, x+y>=3 -> x=3


def fx(kind="", f1="", f2="", f3="", f4="", f5=""):
    """One fixed-column data line: fields at columns 2-3, 5-12, 15-22, 25-36, 40-47, 50-61 (1-based)."""
    line = " " + kind.ljust(2) + " " + f1.ljust(8) + "  " + f2.ljust(8) + "  " + f3.rjust(12)
    if f4:
        line += "   " + f4.ljust(8) + "  " + f5.rjust(12)
    return line.rstrip()


CASES = []


def case(cid, group, what, text, spec, expect=None, why=""):
    CASES.append(dict(id=cid, group=group, what=what, text=text, spec=spec, expect=expect, why=why))


def mutate(old, new, text=BASE):
    assert old in text, old
    return text.replace(old, new)


# ---------------------------------------------------------------------------------------------- format
case("free_control", "format", "control: plain free format", BASE, "accept", OPT3, "baseline")
fixed_text = "\n".join(["NAME          T", "ROWS", fx("N", "obj"), fx("G", "r1"), "COLUMNS",
                        fx("", "x", "obj", "1", "r1", "1"), fx("", "y", "obj", "2", "r1", "1"), "RHS",
                        fx("", "rhs", "r1", "3"), "ENDATA"]) + "\n"
case("fixed_aligned", "format", "strict fixed columns, names <= 8 chars", fixed_text, "accept", OPT3,
     "classic fixed format; whitespace-splitting readers also read it")
sp_text = "\n".join(["NAME          T", "ROWS", fx("N", "OBJ ROW"), fx("G", "MY ROW"), "COLUMNS",
                     fx("", "COL ONE", "OBJ ROW", "1", "MY ROW", "1"), fx("", "COL TWO", "OBJ ROW", "2", "MY ROW", "1"),
                     "RHS", fx("", "RHS", "MY ROW", "3"), "ENDATA"]) + "\n"
case("fixed_names_with_spaces", "format", "fixed columns, names containing spaces", sp_text, "accept", OPT3,
     "spaces inside names are legal in fixed format (field positions, not whitespace, delimit fields)")
case("fixed_shifted_fields", "format", "fixed-format file whose fields are shifted one column but whitespace separated",
     mutate(" x obj 1 r1 1", "  x   obj   1   r1   1"), "either", OPT3,
     "strictly the fields are out of position; whitespace-tolerant readers read it")
case("rhs_without_set_name", "format", "RHS line with the set name omitted (row value)", mutate(" rhs r1 3", " r1 3"), "accept", OPT3,
     "the RHS set name is optional in common practice (and in the fixed-format blank field)")
case("rhs_two_pairs_one_line", "format", "RHS with two row/value pairs on one line",
     mutate("RHS\n rhs r1 3", "RHS\n rhs r1 3 r2 0").replace(" G r1\n", " G r1\n G r2\n").replace("x obj 1 r1 1", "x obj 1 r1 1\n x r2 1"),
     "accept", OPT3, "two pairs per RHS line are standard")
case("lower_case_row_type", "format", "row type in lower case ('g')", mutate(" G r1", " g r1"), "either", OPT3,
     "format says upper case; readers differ")

# ------------------------------------------------------------------------------------------- structure
case("missing_endata", "structure", "no ENDATA line", mutate("ENDATA\n", ""), "reject", None,
     "ENDATA terminates the file; many readers tolerate EOF, the format requires it (judgement call)")
case("text_after_endata", "structure", "garbage after ENDATA", BASE + "this is not mps\n", "accept", OPT3, "nothing after ENDATA is read")
case("missing_name_line", "structure", "no NAME line", mutate("NAME T\n", ""), "either", OPT3, "NAME is documentary; readers differ")
case("empty_file", "structure", "empty file", "", "reject", None, "no model")
case("whitespace_only_file", "structure", "only blank lines", "\n   \n\n", "reject", None, "no model")
case("no_columns_section", "structure", "no COLUMNS section", "NAME T\nROWS\n N obj\nENDATA\n", "either", None,
     "a model with no columns is degenerate; readers differ (accepting it gives an empty model)")
case("no_rhs_section", "structure", "no RHS section (all zero)", mutate("RHS\n rhs r1 3\n", ""), "accept", ("optimal", 0.0), "RHS is optional")
case("bounds_before_columns", "structure", "BOUNDS section before COLUMNS",
     BASE.replace("COLUMNS", "BOUNDS\n UP bnd x 1\nCOLUMNS"), "reject", None, "section order is fixed")
case("unknown_section", "structure", "unknown section name", mutate("RHS", "FOO\n bar baz 1\nRHS"), "reject", None, "unknown section")
case("no_objective_row", "structure", "no N row at all", BASE.replace(" N obj\n", "").replace("x obj 1 r1 1", "x r1 1").replace("y obj 2 r1 1", "y r1 1"),
     "either", ("optimal", 0.0), "a model without an objective row: readers differ (zero objective vs error)")
case("second_objective_row", "structure", "two N rows: the first is the objective, the second is a free row",
     BASE.replace(" G r1", " N free2\n G r1").replace("y obj 2 r1 1", "y obj 2 r1 1\n y free2 7"), "accept", OPT3,
     "extra N rows are discarded by common readers")
case("rows_declared_unused", "structure", "declared row with no entries", BASE.replace(" G r1", " G r1\n L r9"), "accept", OPT3, "an empty row is legal")

# ------------------------------------------------------------------------------------------- duplicates
case("dup_row_name", "duplicates", "two rows with the same name", BASE.replace(" G r1", " G r1\n G r1"), "reject", None,
     "row names are keys")
case("dup_column_noncontiguous", "duplicates", "column x listed, then y, then x again",
     mutate(" y obj 2 r1 1", " y obj 2 r1 1\n x r1 1"), "reject", None,
     "COLUMNS entries of one column must be contiguous (judgement call: some readers merge)")
case("dup_entry_same_cell", "duplicates", "the same (row, column) listed twice", mutate(" x obj 1 r1 1", " x obj 1 r1 1\n x r1 1"),
     "either", None, "format silent: readers sum, keep last, or error")
case("dup_rhs_entry", "duplicates", "RHS for the same row twice", mutate(" rhs r1 3", " rhs r1 3\n rhs r1 5"), "either", None,
     "format silent: last wins or error")
case("dup_bound_same_type", "duplicates", "UP bound listed twice on one column",
     BASE.replace("ENDATA", "BOUNDS\n UP bnd x 10\n UP bnd x 2\nENDATA"), "either", None, "last wins is common; format silent")
case("dup_section_columns", "duplicates", "COLUMNS section twice", mutate("RHS", "COLUMNS\n z obj 1 r1 1\nRHS"), "reject", None,
     "sections occur once")
case("row_and_column_same_name", "duplicates", "a row and a column both named x", BASE.replace(" G r1", " G x"
     ).replace("r1", "x").replace(" x obj 1 x 1", " x obj 1 x 1"), "accept", OPT3, "rows and columns are separate namespaces")
case("name_equals_section_keyword", "duplicates", "a column named RHS and a row named BOUNDS", BASE.replace("r1", "BOUNDS").replace(" y obj", " RHS obj"
     ).replace("rhs BOUNDS 3", "rhs BOUNDS 3"), "accept", ("optimal", 3.0), "data lines are indented, so keywords are not confused with names")

# ------------------------------------------------------------------------------------------ whitespace
case("tabs_free", "whitespace", "tab-separated fields", BASE.replace(" ", "\t"), "either", OPT3, "tabs are not defined by the format; readers differ")
case("tab_leading", "whitespace", "data lines start with a tab", BASE.replace("\n x", "\n\tx").replace("\n y", "\n\ty").replace("\n r", "\n\tr"
     ), "either", OPT3, "indentation by tab")
case("crlf", "whitespace", "CRLF line endings", BASE.replace("\n", "\r\n"), "accept", OPT3, "files written on Windows")
case("cr_only", "whitespace", "CR-only line endings", BASE.replace("\n", "\r"), "either", OPT3, "old Mac files; rare")
case("no_final_newline", "whitespace", "no newline after ENDATA", BASE.rstrip("\n"), "accept", OPT3, "EOF after ENDATA")
case("trailing_spaces", "whitespace", "trailing blanks on every line", BASE.replace("\n", "   \n"), "accept", OPT3, "blanks are insignificant")
case("utf8_bom", "whitespace", "UTF-8 byte-order mark before NAME", b"\xef\xbb\xbf" + BASE.encode(), "either", OPT3,
     "BOM is not part of the format; readers differ")
case("blank_lines_inside", "whitespace", "blank lines between and inside sections", BASE.replace("COLUMNS\n", "COLUMNS\n\n").replace("RHS\n", "\nRHS\n\n"),
     "accept", OPT3, "blank lines are ignored")
case("full_line_comments", "whitespace", "'*' comment lines everywhere", "* head\n" + BASE.replace("ROWS\n", "ROWS\n* c1\n").replace("COLUMNS\n", "* c2\nCOLUMNS\n* c3\n"),
     "accept", OPT3, "a * in column 1 starts a comment")
case("indented_star_comment", "whitespace", "comment with leading space ' * text' inside COLUMNS", mutate("COLUMNS\n", "COLUMNS\n * a comment\n"),
     "reject", None, "comments need the * in column 1; an indented * line is a malformed data line (judgement call)")
case("dollar_inline_comment", "whitespace", "'$' inline comment after the data", mutate(" x obj 1 r1 1", " x obj 1 r1 1 $ note"), "either", OPT3,
     "non-standard extension in some readers")
case("long_names_40", "whitespace", "40-character names in free format", BASE.replace("r1", "r" * 40).replace(" x ", " " + "x" * 40 + " "), "accept", OPT3,
     "free format has no 8-character limit")
case("long_names_300", "whitespace", "300-character names in free format", BASE.replace("r1", "r" * 300), "either", OPT3,
     "no limit in the format; some readers have fixed buffers")
case("names_special_chars", "whitespace", "names such as x[1], a.b-c, 7up", BASE.replace(" x ", " x[1] ").replace(" y ", " 7up ").replace("r1", "a.b-c"),
     "accept", OPT3, "any non-blank characters")
case("very_long_line", "whitespace", "one data line with 100000 spaces of padding", mutate(" x obj 1 r1 1", " x obj 1 r1" + " " * 100000 + "1"), "accept", OPT3,
     "padding is insignificant")

# --------------------------------------------------------------------------------------------- numbers
NUMEXP = {"fortran_D_exponent": 3.0, "plus_sign": 3.0, "leading_dot": 6.0, "trailing_dot": 3.0}
for tag, num_, spec, why in [("fortran_D_exponent", "1D0", "either", "Fortran exponent letter is accepted by several readers"),
                             ("plus_sign", "+1", "accept", "explicit sign"), ("leading_dot", ".5E1", "accept", "x costs 5, y costs 2: y=3, obj 6"),
                             ("trailing_dot", "1.", "accept", "1.0"), ("empty_exponent", "1e", "reject", "malformed"),
                             ("nan_value", "nan", "reject", "not a number in the format"),
                             ("hex_value", "0x10", "reject", "not decimal"), ("decimal_comma", "1,5", "reject", "malformed"),
                             ("trailing_garbage", "3abc", "reject", "malformed"), ("double_sign", "--3", "reject", "malformed"),
                             ("inf_value", "inf", "either", "spelling of infinity as a coefficient is not defined"),
                             ("overflow_1e400", "1e400", "either", "overflows double; readers differ")]:
    case("num_" + tag, "numbers", "coefficient written as '%s'" % num_, mutate(" x obj 1 r1 1", " x obj %s r1 1" % num_), spec,
         ("optimal", NUMEXP[tag]) if tag in NUMEXP else None, why)

# --------------------------------------------------------------------------------------------- bounds
case("bound_undeclared_column", "bounds", "BOUNDS line for a column that is not in COLUMNS (known case)",
     BASE.replace("ENDATA", "BOUNDS\n UP bnd ghost 5\nENDATA"), "reject", None,
     "a bound refers to a column; none was declared (CPLEX and Gurobi report an error; HiGHS ignores the line). Judgement call, not a hard format rule")
case("bound_unknown_type", "bounds", "bound type XX", BASE.replace("ENDATA", "BOUNDS\n XX bnd x 5\nENDATA"), "reject", None, "unknown type")
case("bound_lower_case_type", "bounds", "bound type 'up'", BASE.replace("ENDATA", "BOUNDS\n up bnd x 5\nENDATA"), "either", OPT3, "case rules differ")
case("bound_up_negative_default_lo", "bounds", "UP -1 on a column with default lower bound 0",
     BASE.replace("ENDATA", "BOUNDS\n UP bnd x -1\nENDATA"), "either", None,
     "the format is ambiguous: lower bound -inf with a warning (CPLEX) vs infeasible (lower stays 0) vs error")
case("bound_mi_default_upper", "bounds", "MI on a column (upper bound afterwards?)", BASE.replace("ENDATA", "BOUNDS\n MI bnd x\nENDATA"), "either", None,
     "original IBM reading: upper bound becomes 0; modern readers keep it +inf")
case("bound_bv_with_value", "bounds", "BV with an (ignored) value", BASE.replace("ENDATA", "BOUNDS\n BV bnd x 1\nENDATA"), "accept", ("optimal", 5.0),
     "x binary, y continuous: x=1, y=2 -> obj 5")
case("bound_lo_gt_up", "bounds", "LO 5 and UP 3 on one column", BASE.replace("ENDATA", "BOUNDS\n LO bnd x 5\n UP bnd x 3\nENDATA"), "accept", ("infeasible", None),
     "an empty bound interval is a legal (infeasible) model")
case("bound_two_sets", "bounds", "two named BOUNDS sets (second sets x <= 1)", BASE.replace("ENDATA", "BOUNDS\n UP b1 x 10\n UP b2 x 1\nENDATA"), "either", None,
     "only one set is selected by the user; readers use the first or apply all")
case("bound_sc_semicontinuous", "bounds", "SC (semi-continuous) bound", BASE.replace("ENDATA", "BOUNDS\n SC bnd x 10\nENDATA"), "either", None,
     "extension, not in the basic format")
case("bound_free_then_up", "bounds", "FR followed by UP 2 on x", BASE.replace("ENDATA", "BOUNDS\n FR bnd x\n UP bnd x 2\nENDATA"), "accept", ("optimal", 4.0),
     "sequential application: x in (-inf, 2], y >= 0: min x+2y s.t. x+y>=3 -> x=2, y=1, obj 4")
MK = "    MARKER                 'MARKER'                 '%s'\n"
case("integer_marker_no_bounds", "bounds", "integer column in a MARKER block with no bound, minimising -x",
     "NAME T\nROWS\n N obj\n L r1\nCOLUMNS\n" + MK % "INTORG" + " x obj -1 r1 0\n" + MK % "INTEND" + "RHS\nENDATA\n",
     "either", None, "default upper bound of an integer column: +inf (modern) or 1 (historic): unbounded vs optimal -1")
case("marker_unbalanced", "bounds", "INTORG marker never closed", BASE.replace("COLUMNS\n", "COLUMNS\n" + MK % "INTORG"), "either", OPT3,
     "an unclosed block is malformed in principle; readers treat the rest as integer")
case("marker_bad_token", "bounds", "MARKER line with an unknown type", BASE.replace("COLUMNS\n", "COLUMNS\n" + MK % "INTFOO"), "reject", None, "malformed marker")

case("bound_missing_value", "bounds", "UP bound with no value", BASE.replace("ENDATA", "BOUNDS\n UP bnd x\nENDATA"), "reject", None, "UP needs a value")
case("bound_bad_number", "bounds", "UP bound with a non-number", BASE.replace("ENDATA", "BOUNDS\n UP bnd x abc\nENDATA"), "reject", None, "malformed number")

# ---------------------------------------------------------------------------------------- rhs/ranges
case("rhs_missing_value", "rhs_ranges", "RHS entry with a row but no value", mutate(" rhs r1 3", " rhs r1"), "reject", None, "missing value")
case("ranges_unknown_row", "rhs_ranges", "RANGES for a row that does not exist", mutate("ENDATA", "RANGES\n rng nosuch 4\nENDATA"), "reject", None, "unknown row")
case("rows_line_three_fields", "rhs_ranges", "ROWS line with a third field", mutate(" G r1", " G r1 extra"), "reject", None, "malformed row declaration (fixed format would read the extra text as nothing)")
case("empty_optional_sections", "rhs_ranges", "empty RHS, RANGES and BOUNDS sections", mutate("RHS\n rhs r1 3\n", "RHS\nRANGES\nBOUNDS\n"), "accept", ("optimal", 0.0), "empty sections are legal")
case("name_line_extra_words", "rhs_ranges", "NAME line with several words", mutate("NAME T", "NAME  my model name"), "accept", OPT3, "NAME is documentary")
case("rhs_unknown_row", "rhs_ranges", "RHS for a row that does not exist", mutate(" rhs r1 3", " rhs r1 3\n rhs nosuch 4"), "reject", None, "unknown row")
case("columns_unknown_row", "rhs_ranges", "COLUMNS entry for an undeclared row", mutate(" x obj 1 r1 1", " x obj 1 r1 1 nosuch 2"), "reject", None, "unknown row")
case("range_negative_on_E", "rhs_ranges", "E row with R<0: interval [rhs+R, rhs]",
     "NAME T\nROWS\n N obj\n E r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 5\nRANGES\n rng r1 -2\nENDATA\n", "accept", ("optimal", 3.0), "x in [3,5], min x = 3")
case("range_on_G_negative", "rhs_ranges", "G row with R<0: |R| is used, interval [rhs, rhs+|R|]",
     "NAME T\nROWS\n N obj\n G r1\nCOLUMNS\n x obj -1 r1 1\nRHS\n rhs r1 5\nRANGES\n rng r1 -2\nENDATA\n", "accept", ("optimal", -7.0), "x in [5,7], min -x = -7")
case("range_on_N_row", "rhs_ranges", "RANGES entry for the objective row (ignored)", mutate("ENDATA", "RANGES\n rng obj 4\nENDATA"), "accept", OPT3, "ranges on N rows have no meaning")
case("range_zero", "rhs_ranges", "R = 0 on an L row (becomes equality)",
     "NAME T\nROWS\n N obj\n L r1\nCOLUMNS\n x obj -1 r1 1\nRHS\n rhs r1 5\nRANGES\n rng r1 0\nENDATA\n", "accept", ("optimal", -5.0), "x in [5,5]")
case("rhs_on_objective_row", "rhs_ranges", "RHS on the objective row (constant term, sign: obj = c'x - rhs)", mutate(" rhs r1 3", " rhs r1 3\n rhs obj 10"), "accept", ("optimal", -7.0),
     "objective constant is minus the RHS entry (standard; HiGHS, CPLEX, Gurobi agree)")
case("objsense_max_section", "rhs_ranges", "OBJSENSE section with MAX on its own line",
     "NAME T\nOBJSENSE\n    MAX\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n", "accept", ("optimal", 4.0), "maximise x s.t. x<=4")
case("objsense_inline", "rhs_ranges", "'OBJSENSE MAX' on one line",
     "NAME T\nOBJSENSE MAX\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n", "either", ("optimal", 4.0), "inline form is common but not universal")
case("objsense_maximize", "rhs_ranges", "OBJSENSE with the word MAXIMIZE",
     "NAME T\nOBJSENSE\n    MAXIMIZE\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n", "accept", ("optimal", 4.0), "MAXIMIZE is the long form")
case("objsense_lowercase", "rhs_ranges", "OBJSENSE with lower-case 'max'",
     "NAME T\nOBJSENSE\n    max\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n", "either", ("optimal", 4.0), "case rules differ")
case("objsense_garbage", "rhs_ranges", "OBJSENSE with an unknown word",
     "NAME T\nOBJSENSE\n    SIDEWAYS\nROWS\n N obj\n L r1\nCOLUMNS\n x obj 1 r1 1\nRHS\n rhs r1 4\nENDATA\n", "reject", None, "unknown sense")
case("columns_line_four_tokens", "rhs_ranges", "COLUMNS line with an odd number of fields", mutate(" x obj 1 r1 1", " x obj 1 r1"), "reject", None, "missing value")
case("columns_line_seven_tokens", "rhs_ranges", "COLUMNS line with three pairs", mutate(" x obj 1 r1 1", " x obj 1 r1 1 obj 1"), "reject", None, "at most two pairs per line")


def run_engine(engine, path):
    js = path + ".json"
    try:
        p = subprocess.run([engine, path, "--time-limit", "10", "--json", js], capture_output=True, timeout=30)
    except Exception as e:
        return dict(read="crash", status="crash", obj=None, msg=repr(e)[:100])
    try:
        j = json.load(open(js))
    except Exception:  # the engine writes bare nan/inf into its JSON for a non-finite objective: invalid JSON
        out = p.stdout.decode(errors="replace").split()
        st = out[1] if len(out) > 1 and out[0] == "status" else "no_status"
        try:
            ob = float(out[3]) if len(out) > 3 and out[2] == "objective" else None
        except ValueError:
            ob = None
        return dict(read="accept", status=st, obj=ob, msg="invalid JSON written (rc=%s)" % p.returncode, invalid_json=True)
    st = j.get("status")
    return dict(read="reject" if st == "parse_error" else "accept", status=st, obj=j.get("objective"), msg=(j.get("message") or "")[:100])


def run_highs(path):
    try:
        h, ok, rs = read_ref(path)
    except Exception as e:
        return dict(read="reject", status="read_exception", obj=None, msg=repr(e)[:100])
    if not ok:
        return dict(read="reject", status="read_error", obj=None, msg=rs)
    try:
        import highspy
        h.setOptionValue("presolve", "off")
        h.setOptionValue("time_limit", 10.0)
        h.run()
        ms = h.getModelStatus()
        from ref import STATUS
        st = STATUS.get(ms, str(ms).split(".")[-1])
        obj = h.getInfo().objective_function_value if st == "optimal" else None
    except Exception as e:
        return dict(read="accept", status="solve_exception", obj=None, msg=repr(e)[:100])
    return dict(read="accept", status=st, obj=obj, msg=rs)


def matches(res, expect):
    if expect is None:
        return True
    st, obj = expect
    if res["status"] != st:
        return False
    return obj is None or (res["obj"] is not None and abs(res["obj"] - obj) <= 1e-6 * max(1, abs(obj)))


def describe(r):
    if r["read"] == "reject":
        return "error"
    if r["status"] == "optimal" and r["obj"] is not None:
        return "optimal %.6g" % r["obj"]  # nan prints as 'nan'

    return r["status"]


def verdict(c, o, h):
    spec = c["spec"]
    if o["status"] in ("crash", "no_status") or (o["read"] == "accept" and o["status"] == "optimal" and o["obj"] is not None
                                                  and o["obj"] != o["obj"]) or o.get("invalid_json"):
        return "FAIL", "n/a (engine reports a non-finite or unreadable result whatever the input convention)"
    if spec == "accept":
        oo, hh = o["read"] == "accept" and matches(o, c["expect"]), h["read"] == "accept" and matches(h, c["expect"])
        who = "both" if oo and hh else "ours only" if oo else "HiGHS only" if hh else "neither"
        return ("PASS" if oo else "FAIL"), who
    if spec == "reject":
        oo, hh = o["read"] == "reject", h["read"] == "reject"
        who = "both" if oo and hh else "ours only" if oo else "HiGHS only" if hh else "neither"
        return ("PASS" if oo else "FAIL"), who
    same = describe(o) == describe(h)
    return "PASS", ("same behaviour" if same else "spec silent: both defensible")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--known", default="tests/known_failures.tsv")
    a = ap.parse_args()
    d = os.path.join(a.out, "parser")
    os.makedirs(d, exist_ok=True)
    rows = []
    for c in CASES:
        p = os.path.join(d, c["id"] + ".mps")
        data = c["text"] if isinstance(c["text"], bytes) else c["text"].encode()
        with open(p, "wb") as f:
            f.write(data)
        o, h = run_engine(os.path.abspath(a.engine), p), run_highs(p)
        v, who = verdict(c, o, h)
        rows.append(dict(id=c["id"], group=c["group"], what=c["what"], spec=c["spec"], ours=describe(o), highs=describe(h),
                         verdict=v, right=who, why=c["why"], ours_msg=o["msg"]))
    json.dump(rows, open(os.path.join(a.out, "parser_table.json"), "w"), indent=1)
    with open(os.path.join(a.out, "parser_table.md"), "w") as f:
        f.write("| input | what | spec says | ours | HiGHS | who is right | ours verdict |\n|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write("| `%s` | %s | %s | %s | %s | %s | %s |\n" % (r["id"], r["what"].replace("|", "/"), r["spec"], r["ours"], r["highs"], r["right"], r["verdict"]))
    n = len(rows)
    npass = sum(r["verdict"] == "PASS" for r in rows)
    known = set()
    kf = a.known
    if os.path.exists(kf):
        for ln in open(kf):
            if ln.startswith("parser:"):
                known.add(ln.split("\t")[0][len("parser:"):])
    fails = [r for r in rows if r["verdict"] == "FAIL"]
    new = [r for r in fails if r["id"] not in known]
    print("parser table: %d inputs, %d pass, %d fail (%d known, %d new); spec-silent inputs count as pass" % (
        n, npass, len(fails), len(fails) - len(new), len(new)))
    for r in fails:
        print("  %-6s %-30s ours=%-14s HiGHS=%-14s" % ("KNOWN" if r["id"] in known else "NEW", r["id"], r["ours"], r["highs"]))
    if new:
        sys.exit(1)


if __name__ == "__main__":
    main()
