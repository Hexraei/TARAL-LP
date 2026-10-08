#!/usr/bin/env python3
"""Tests for rim_vector_report.py on generated MPS files (no solver involved)."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rim_vector_report as rr

TMP = tempfile.mkdtemp()
fails = []

def expect(name, cond, got=None):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else " got=" + repr(got)))
    if not cond:
        fails.append(name)

def write(name, text):
    p = os.path.join(TMP, name + ".mps")
    with open(p, "w") as f:
        f.write(text)
    return p

BASE = """NAME          T
ROWS
 N  COST
 G  R1
 L  R2
COLUMNS
    X         COST      1.0       R1        1.0
    X         R2        1.0
"""
BOUNDS1 = """BOUNDS
 UP BND       X         10.
"""
END = "ENDATA\n"

# 1. single named set: selected, nothing skipped
r = rr.scan(write("single", BASE + "RHS\n    B         R1        5.0\n" + BOUNDS1 + END))
occ = r["rhs"]["occurrences"][0]
expect("single_set_selected", occ["selected"] == "B" and occ["name_skipped"] == 0
       and occ["entries_by_name"] == {"B": 1}, occ)
expect("single_set_bounds_selected", r["bounds"]["occurrences"][0]["selected"] == "BND")
expect("single_set_ranges_absent", r["ranges"] is None)

# 2. interleaved foreign name first in one section (supplied-reproducer shape)
r = rr.scan(write("foreign", BASE + "RHS\n    RHS       COST      -100.0\n    B         R1        5.0\n" + END))
occ = r["rhs"]["occurrences"][0]
expect("foreign_selects_first_name", occ["selected"] == "RHS", occ)
expect("foreign_skips_other_name", occ["name_skipped"] == 1
       and occ["entries_by_name"] == {"RHS": 1, "B": 1}, occ)

# 3. two-pair data lines: pairs counted individually
r = rr.scan(write("twopair", BASE + "RHS\n    B         R1        5.0       R2        7.0\n"
                  "    F         R2        9.0\n" + END))
occ = r["rhs"]["occurrences"][0]
expect("twopair_counts_pairs", occ["entries_by_name"] == {"B": 2, "F": 1}
       and occ["name_skipped"] == 1, occ)

# 4. nameless lines (even token count) counted separately, never skipped
r = rr.scan(write("nameless", BASE + "RHS\n    B         R1        5.0\n"
                  "    R2        3.0\n" + END))
occ = r["rhs"]["occurrences"][0]
expect("nameless_counted_applied", occ["nameless"] == 1 and occ["name_skipped"] == 0
       and occ["selected"] == "B", occ)

# 5. three names interleaved: both later names skipped
r = rr.scan(write("threenames", BASE + "RHS\n    A         R1        1.0\n"
                  "    B         R1        2.0\n    C         R2        3.0\n    A         R2        4.0\n" + END))
occ = r["rhs"]["occurrences"][0]
expect("three_names_skip_two", occ["selected"] == "A" and occ["name_skipped"] == 2
       and occ["entries_by_name"] == {"A": 2, "B": 1, "C": 1}, occ)

# 6. RANGES selection is independent per kind
r = rr.scan(write("ranges", BASE + "RHS\n    B         R1        5.0\n"
                  "RANGES\n    RNG1      R1        2.0\n    RNG2      R2        3.0\n" + BOUNDS1 + END))
expect("ranges_independent_selection", r["ranges"]["occurrences"][0]["selected"] == "RNG1"
       and r["ranges"]["occurrences"][0]["name_skipped"] == 1
       and r["rhs"]["occurrences"][0]["selected"] == "B", r["ranges"])

# 7. BOUNDS: two names, no-value types, and 3-token nameless line
r = rr.scan(write("bounds", BASE + "RHS\n    B         R1        5.0\n"
                  "BOUNDS\n UP BND       X         10.\n FX ALT       X         4.\n"
                  " MI BND       X\n UP X         8.\n" + END))
occ = r["bounds"]["occurrences"][0]
expect("bounds_selection_and_skip", occ["selected"] == "BND" and occ["name_skipped"] == 1
       and occ["entries_by_name"] == {"BND": 2, "ALT": 1}, occ)
expect("bounds_nameless_lines", occ["nameless"] == 1, occ)  # "UP X 8." has no set field

# 8. repeated RHS section header: noncanonical, semantics not assumed
r = rr.scan(write("repeated_rhs", BASE + "RHS\n    B         R1        5.0\n"
                  "RHS\n    C         R2        1.0\n" + END))
expect("repeated_rhs_labelled", r["rhs"]["repeated"] is True
       and r["rhs"]["semantics"] == "not_assumed"
       and len(r["rhs"]["occurrences"]) == 2, r["rhs"])

# 9. repeated BOUNDS section header: same treatment
r = rr.scan(write("repeated_bounds", BASE + "RHS\n    B         R1        5.0\n"
                  + BOUNDS1 + "BOUNDS\n LO ALT       X         1.\n" + END))
expect("repeated_bounds_labelled", r["bounds"]["repeated"] is True
       and r["bounds"]["semantics"] == "not_assumed"
       and len(r["bounds"]["occurrences"]) == 2, r["bounds"])

# 10. no rim sections at all
r = rr.scan(write("norim", BASE + END))
expect("no_rim_sections_null", r["rhs"] is None and r["ranges"] is None and r["bounds"] is None)

# 11. indented data line starting with a keyword token is not a header (regression)
r = rr.scan(write("keyword_row", """NAME          T
ROWS
 N  COST
 G  RHS
COLUMNS
    X         COST      1.0       RHS       1.0
RHS
    B         RHS       5.0
ENDATA
"""))
occ = r["rhs"]["occurrences"][0]
expect("indented_keyword_token_not_header", occ["selected"] == "B"
       and occ["entries_by_name"] == {"B": 1} and r["rhs"]["repeated"] is False, r["rhs"])

print("FAILED %d" % len(fails) if fails else "ALL PASS")
sys.exit(1 if fails else 0)
