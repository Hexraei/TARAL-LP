#!/usr/bin/env python3
"""Tests for columns_structure_report.py on generated MPS files (no solver involved)."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import columns_structure_report as cr

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

HDR = "NAME          T\nROWS\n N  COST\n G  R1\n L  R2\n"
END = "ENDATA\n"

# 1. clean file: counts exact, no flags, status clean
r = cr.scan(write("clean", HDR + "COLUMNS\n    X         COST      1.0       R1        1.0\n"
                  "    Y         R2        1.0\n" + END))
c = r["columns"]
expect("clean_counts", c["distinct"] == 2 and c["data_lines"] == 2
       and c["entries_total"] == 3 and c["matrix_entries"] == 2
       and c["objective_row_entries"] == 1, c)
expect("clean_status", r["tallies_status"] == "clean"
       and r["repeated_section_headers"] == {"ROWS": 0, "COLUMNS": 0}, r["tallies_status"])

# 2. objective row is the first N row; later N rows are extra free
r = cr.scan(write("two_n", "ROWS\n N  OBJ1\n N  OBJ2\n G  R1\nCOLUMNS\n"
                  "    X         OBJ1      1.0       OBJ2      2.0\n    X         R1        1.0\n" + END))
expect("first_n_objective", r["rows"]["objective_row"] == "OBJ1"
       and r["rows"]["extra_free_rows"] == ["OBJ2"], r["rows"])
expect("free_n_dropped", r["columns"]["free_n_row_entries_dropped"] == 1
       and r["columns"]["objective_row_entries"] == 1
       and r["columns"]["matrix_entries"] == 1, r["columns"])

# 3. duplicate (column, row) entries counted per category: matrix / objective / free-N
r = cr.scan(write("dup", "ROWS\n N  COST\n N  FREE2\n G  R1\nCOLUMNS\n"
                  "    X         COST      1.0       R1        1.0\n"
                  "    X         R1        2.0\n"
                  "    X         COST      3.0       FREE2     1.0\n"
                  "    X         FREE2     2.0\n" + END))
c = r["columns"]
expect("dup_matrix", c["duplicate_matrix_entries"]["count"] == 1
       and c["duplicate_matrix_entries"]["examples"][0]["row"] == "R1",
       c["duplicate_matrix_entries"])
expect("dup_objective", c["duplicate_objective_row_entries"]["count"] == 1
       and c["duplicate_objective_row_entries"]["examples"][0]["row"] == "COST",
       c["duplicate_objective_row_entries"])
expect("dup_free_n", c["duplicate_free_n_row_entries"]["count"] == 1
       and c["duplicate_free_n_row_entries"]["examples"][0]["row"] == "FREE2",
       c["duplicate_free_n_row_entries"])

# 4. same column in two separate blocks
r = cr.scan(write("twoblocks", HDR + "COLUMNS\n    X         R1        1.0\n"
                  "    Y         R1        1.0\n    X         R2        1.0\n" + END))
nb = r["columns"]["nonconsecutive_column_blocks"]
expect("two_blocks_flagged", nb["count"] == 1 and nb["columns"] == ["X"]
       and r["tallies_status"].startswith("hypothetical"), nb)

# 5. bad token counts (2, 4, 6 tokens) counted, not interpreted
r = cr.scan(write("badtok", HDR + "COLUMNS\n    X         R1\n"
                  "    Y         R1        1.0       R2\n"
                  "    Z         R1        1.0       R2        1.0       R2        2.0\n"
                  "    W         R1        1.0\n" + END))
bt = r["columns"]["bad_token_count_lines"]
expect("bad_token_counts", bt["count"] == 3 and r["columns"]["data_lines"] == 1, bt)

# 5b. an unreadable line ends a column's consecutive run
r = cr.scan(write("badtok_block", HDR + "COLUMNS\n    X         R1        1.0\n"
                  "    Y         R1\n    X         R2        1.0\n" + END))
expect("bad_line_breaks_block", r["columns"]["nonconsecutive_column_blocks"]["count"] == 1
       and r["columns"]["nonconsecutive_column_blocks"]["columns"] == ["X"],
       r["columns"]["nonconsecutive_column_blocks"])

# 6. unknown row reference counted
r = cr.scan(write("unkrow", HDR + "COLUMNS\n    X         R9        1.0       R1        1.0\n" + END))
expect("unknown_row", r["columns"]["unknown_row_references"] == 1
       and r["columns"]["matrix_entries"] == 1, r["columns"])

# 7. marker region: columns between INTORG and INTEND marked integer
r = cr.scan(write("markers", HDR + "COLUMNS\n    X         R1        1.0\n"
                  "    MARKER                 'MARKER'                 'INTORG'\n"
                  "    Y         R1        1.0\n    Z         R2        1.0\n"
                  "    MARKER                 'MARKER'                 'INTEND'\n"
                  "    W         R1        1.0\n" + END))
m = r["columns"]["markers"]
expect("marker_region", m["intorg"] == 1 and m["intend"] == 1
       and m["columns_marked_integer"] == 2
       and m["open_at_columns_section_end"] is False, m)

# 8. INTORG open at section end: open-ended columns are only the final region's
r = cr.scan(write("openint", HDR + "COLUMNS\n"
                  "    MARKER                 'MARKER'                 'INTORG'\n"
                  "    Y         R1        1.0\n"
                  "    MARKER                 'MARKER'                 'INTEND'\n"
                  "    W         R1        1.0\n"
                  "    MARKER                 'MARKER'                 'INTORG'\n"
                  "    Q         R1        1.0\n    P         R2        1.0\n"
                  "RHS\n    B         R1        5.0\n" + END))
m = r["columns"]["markers"]
expect("open_at_end", m["open_at_columns_section_end"] is True
       and m["columns_marked_integer"] == 3
       and m["columns_marked_integer_open_ended"] == 2, m)

# 9. ambiguous marker lines split: both-present vs neither-present; INTEND without open
r = cr.scan(write("badmarker", HDR + "COLUMNS\n"
                  "    MARKER                 'MARKER'    'INTORG'    'INTEND'\n"
                  "    MARKER                 'MARKER'                 'INTFLD'\n"
                  "    MARKER                 'MARKER'                 'INTEND'\n"
                  "    X         R1        1.0\n" + END))
m = r["columns"]["markers"]
expect("bad_markers", m["marker_line_intorg_and_intend"] == 1
       and m["marker_line_no_intorg_intend"] == 1
       and m["intend_without_open"] == 1
       and m["columns_marked_integer"] == 0, m)

# 9b. redundant INTORG while open: counted, region preserved (not restarted)
r = cr.scan(write("redintorg", HDR + "COLUMNS\n    X         R1        1.0\n"
                  "    MARKER                 'MARKER'                 'INTORG'\n"
                  "    Y         R1        1.0\n"
                  "    MARKER                 'MARKER'                 'INTORG'\n"
                  "    Z         R2        1.0\n"
                  "RHS\n    B         R1        5.0\n" + END))
m = r["columns"]["markers"]
expect("redundant_intorg_region_preserved", m["redundant_intorg_while_open"] == 1
       and m["columns_marked_integer"] == 2
       and m["columns_marked_integer_open_ended"] == 2
       and m["open_at_columns_section_end"] is True, m)

# 10. a column literally named MARKER without quotes is data, not a marker line
r = cr.scan(write("notmarker", HDR + "COLUMNS\n    MARKER     R1        1.0\n" + END))
expect("quoted_marker_only", r["columns"]["markers"]["marker_lines"] == 0
       and r["columns"]["distinct"] == 1, r["columns"]["markers"])

# 11. repeated COLUMNS header counted
r = cr.scan(write("repsec", HDR + "COLUMNS\n    X         R1        1.0\n"
                  "COLUMNS\n    Y         R1        1.0\n" + END))
expect("repeated_columns", r["repeated_section_headers"]["COLUMNS"] == 1
       and r["columns"]["distinct"] == 2, r["repeated_section_headers"])

# 12. ENDATA stops the scan; trailing COLUMNS content ignored
r = cr.scan(write("trailing", HDR + "COLUMNS\n    X         R1        1.0\n" + END
                  + "COLUMNS\n    ZZZ       R1        1.0\n"))
expect("endata_stop", r["columns"]["distinct"] == 1, r["columns"])

# 13. comments and blank lines skipped; unknown column-1 header counted, ends section
r = cr.scan(write("comments", "* comment\n" + HDR + "\nCOLUMNS\n* mid comment\n"
                  "    X         R1        1.0\nFOOBAR\n    Y         R1        1.0\n" + END))
expect("comments_unknown_header", r["columns"]["distinct"] == 1
       and r["unknown_column1_headers"] == 1, r["unknown_column1_headers"])
expect("unknown_header_marks_hypothetical", r["tallies_status"].startswith("hypothetical"),
       r["tallies_status"])

# 13b. not_validated list is reported
r = cr.scan(write("clean2", HDR + "COLUMNS\n    X         R1        1.0\n" + END))
expect("not_validated_reported", "numeric token syntax" in r["not_validated"]
       and r["tallies_status"] == "clean", r.get("not_validated"))

# 14. ROWS anomalies: bad token count, bad type, duplicate name each counted
r = cr.scan(write("badrows", "ROWS\n N  COST\n N  COST EXTRA\n Q  R1\n E  COST\n G  R2\n"
                  "COLUMNS\n    X         R2        1.0\n" + END))
expect("bad_rows_tallies", r["rows"]["bad_line_count"] == 1
       and r["rows"]["bad_type_count"] == 1
       and r["rows"]["duplicate_name_count"] == 1
       and r["rows"]["objective_row"] == "COST", r["rows"])

# 15. every symbol_map citation points at src/mps.cpp with a line reference
bad = [k for k, v in cr.SYMBOL_MAP.items() if "src/mps.cpp:" not in v]
expect("symbol_map_cites_lines", bad == [], bad)

print("FAILED %d" % len(fails) if fails else "ALL PASS")
sys.exit(1 if fails else 0)
