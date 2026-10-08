#!/usr/bin/env python3
"""Static MPS COLUMNS-section structural report.

Usage: columns_structure_report.py file.mps [file2.mps ...]

INFORMATIONAL ONLY. This checker reports flags and counts for the ROWS and
COLUMNS sections of an MPS file. It does NOT predict whether any file is
accepted or rejected, and it takes no position on reader semantics: where a
counted construct is one the reader acts on, the reader's behavior is quoted
from the source (symbol and line, src/mps.cpp at commit
dd217e607caef1e56b3a681cbceafd9a0d81c589) and adjudication is left to the
reader's maintainers.

Scope and divergence notes:
- Free-format tokenization is assumed, mirroring tokens(line, section,
  fixed=false) at src/mps.cpp:68-81 (whitespace split, C-locale). The reader
  additionally retries the whole file with fixed-column tokenization when the
  free reading fails (read_mps, src/mps.cpp:363-373); that fallback is NOT
  modelled here. Files written for the fixed-column reading can tokenize
  differently here, and all counts for such files are local hypotheticals.
- This is a STATIC checker: it runs no solver and modifies nothing. It is not
  a general parser and it diverges from the reader on invalid or noncanonical
  input. Repeated ROWS / COLUMNS section headers, bad token counts, unknown
  row references and nonconsecutive column blocks are counted, not
  interpreted; once any of them is present, downstream tallies are labelled
  "hypothetical" because the reader's own parse state at that point is its
  own business, not this report's.
- Scanning stops at ENDATA (src/mps.cpp:162-164); trailing content is ignored
  and not compared against the reader. Blank lines and lines whose first
  character is '*' are skipped (src/mps.cpp:154). Lines starting in column 1
  whose first token is not a known header are counted as
  unknown_column1_headers, not interpreted.

Counted constructs and the reader behavior noted for each (all citations
src/mps.cpp at dd217e607caef1e56b3a681cbceafd9a0d81c589):
- duplicate (column, row) entries inside one consecutive block: summed in
  file order by flush_column, zero sums dropped (flush_column, lines 111-123;
  sum at line 118, zero-drop at line 121).
- entries naming extra free N rows (every N row after the first): dropped;
  row_of returns -2 for them and the entry is neither cost nor matrix
  (row_of, lines 139-144; use at lines 237-241). First N row is the objective
  (line 200); later N rows go to free_rows (line 201; comment line 102).
- entries naming the objective row: added to md.cost (line 240).
- a column listed in two separate (non-consecutive) blocks: ParseError at
  line 233 (lookup path lines 219-234).
- COLUMNS data lines with a token count other than 3 or 5: ParseError at
  line 217.
- MARKER lines (token[1] == "'MARKER'", line 210): 'INTORG' and 'INTEND' on
  the same line is a ParseError (lines 211-213); otherwise in_int is set from
  'INTORG' presence (line 214) and every column declared while in_int holds
  is flagged integer (md.is_int push, line 229). An 'INTORG' still open at
  the end of the COLUMNS section is not an error in the reader; the columns
  it covers are counted here as columns_marked_integer_open_ended.
- repeated ROWS or COLUMNS section headers: ParseError at lines 182-183.
- entries naming rows never declared in ROWS: ParseError from row_of at
  line 144.

Local conventions where the reader stops first (its ParseError ends the
parse, so any state after such a line is this report's own choice, not a
reader behavior): a skipped bad-token-count line ends the column's
consecutive run for block accounting, and scanning continues past
noncanonical constructs so their full counts are visible.
"""
import hashlib, json, sys

KNOWN_HEADERS = {"NAME", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS",
                 "OBJSENSE", "OBJNAME", "QUADOBJ", "QMATRIX", "QSECTION",
                 "QCMATRIX", "SOS", "ENDATA"}
ROW_TYPES = {"N", "E", "L", "G"}
EXAMPLE_CAP = 10

SYMBOL_MAP = {
    "tokenization": "src/mps.cpp:68-81 tokens(line, section, fixed=false)",
    "fixed_format_fallback_not_modelled": "src/mps.cpp:363-373 read_mps",
    "comment_and_blank_skip": "src/mps.cpp:154",
    "endata_stop": "src/mps.cpp:162-164",
    "rows_first_n_is_objective": "src/mps.cpp:200 obj_row",
    "rows_later_n_are_free": "src/mps.cpp:201 free_rows (comment src/mps.cpp:102)",
    "rows_bad_line": "src/mps.cpp:195",
    "rows_duplicate_name": "src/mps.cpp:197-198",
    "rows_bad_type": "src/mps.cpp:207",
    "repeated_rows_or_columns_section": "src/mps.cpp:182-183 seen_sections",
    "marker_line_detection": "src/mps.cpp:210",
    "marker_same_line_intorg_intend": "src/mps.cpp:211-213",
    "marker_sets_in_int": "src/mps.cpp:214 in_int",
    "columns_bad_token_count": "src/mps.cpp:217",
    "column_declared_while_in_int": "src/mps.cpp:229 md.is_int",
    "column_two_separate_blocks": "src/mps.cpp:233",
    "objective_row_entry_to_cost": "src/mps.cpp:240 md.cost",
    "matrix_entry_collected": "src/mps.cpp:241 col_entries",
    "free_n_row_entry_dropped": "src/mps.cpp:139-144 row_of (-2), 237-241",
    "unknown_row_reference": "src/mps.cpp:144 row_of",
    "duplicate_entries_summed": "src/mps.cpp:111-123 flush_column (sum line 118, zero-drop line 121)",
}


def scan(path):
    with open(path, "rb") as f:
        raw = f.read()
    rep = {
        "file": path,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "format_assumed": "free (fixed-column fallback not modelled)",
        "adjudication": "none: flags/counts only; reader behavior quoted from source, not predicted",
        "rows": {"declared_constraint_rows": 0, "objective_row": None,
                 "extra_free_rows": [], "bad_line_count": 0,
                 "bad_type_count": 0, "duplicate_name_count": 0},
        "columns": {"distinct": 0, "data_lines": 0, "entries_total": 0,
                    "matrix_entries": 0, "objective_row_entries": 0,
                    "free_n_row_entries_dropped": 0,
                    "unknown_row_references": 0,
                    "duplicate_pair_entries": {"count": 0, "examples": []},
                    "nonconsecutive_column_blocks": {"count": 0, "columns": []},
                    "bad_token_count_lines": {"count": 0, "lines": []},
                    "markers": {"marker_lines": 0, "intorg": 0, "intend": 0,
                                "same_line_intorg_intend": 0,
                                "intend_without_open": 0,
                                "open_at_columns_section_end": False,
                                "columns_marked_integer": 0,
                                "columns_marked_integer_open_ended": 0}},
        "repeated_section_headers": {"ROWS": 0, "COLUMNS": 0},
        "unknown_column1_headers": 0,
        "tallies_status": "clean",
        "symbol_map": SYMBOL_MAP,
    }
    section = None
    lineno = 0
    seen_sections = set()
    row_kind = {}          # name -> "N" (objective) / "n" (extra free) / constraint type
    obj_row = None
    col_order = []         # distinct column names in first-appearance order
    col_blocks = {}        # column name -> number of separate blocks
    cur_col = None
    cur_block_rows = set()
    in_int = False
    int_columns = set()
    region_columns = []   # columns declared in the current open marker region
    open_ended_int_columns = set()
    noncanonical = False

    def close_block():
        nonlocal cur_col, cur_block_rows
        cur_col, cur_block_rows = None, set()

    for line in raw.decode("utf-8", "replace").splitlines():
        lineno += 1
        if not line.strip() or line.startswith("*"):
            continue
        if line[0] not in " \t":
            head = line.split()[0].upper()
            if head == "ENDATA":
                break
            if head in KNOWN_HEADERS:
                if head in ("ROWS", "COLUMNS"):
                    if head in seen_sections:
                        rep["repeated_section_headers"][head] += 1
                        noncanonical = True
                    seen_sections.add(head)
                if section == "COLUMNS" and head != "COLUMNS":
                    close_block()
                    if in_int:
                        rep["columns"]["markers"]["open_at_columns_section_end"] = True
                        open_ended_int_columns = set(region_columns)
                        region_columns = []
                        in_int = False
                section = head
            else:
                if section == "COLUMNS":
                    close_block()
                    if in_int:
                        rep["columns"]["markers"]["open_at_columns_section_end"] = True
                        open_ended_int_columns = set(region_columns)
                        region_columns = []
                        in_int = False
                section = None
                rep["unknown_column1_headers"] += 1
            continue
        if section == "ROWS":
            t = line.split()
            if len(t) != 2 or not t[1]:
                rep["rows"]["bad_line_count"] += 1
                noncanonical = True
                continue
            typ = t[0][0].upper() if t[0] else "?"
            name = t[1]
            if name in row_kind:
                rep["rows"]["duplicate_name_count"] += 1
                noncanonical = True
                continue
            if typ == "N":
                if obj_row is None:
                    obj_row = name
                    rep["rows"]["objective_row"] = name
                    row_kind[name] = "N"
                else:
                    rep["rows"]["extra_free_rows"].append(name)
                    row_kind[name] = "n"
            elif typ in ("E", "L", "G"):
                row_kind[name] = typ
                rep["rows"]["declared_constraint_rows"] += 1
            else:
                rep["rows"]["bad_type_count"] += 1
                noncanonical = True
            continue
        if section != "COLUMNS":
            continue
        t = line.split()
        c = rep["columns"]
        if len(t) >= 2 and t[1] == "'MARKER'":
            m = c["markers"]
            m["marker_lines"] += 1
            start = "'INTORG'" in t
            end = "'INTEND'" in t
            if start:
                m["intorg"] += 1
            if end:
                m["intend"] += 1
            if start == end:
                m["same_line_intorg_intend"] += 1
                noncanonical = True
            elif end:
                if not in_int:
                    m["intend_without_open"] += 1
                in_int = False
                region_columns = []
            else:
                in_int = True
                region_columns = []
            continue
        if len(t) not in (3, 5):
            c["bad_token_count_lines"]["count"] += 1
            if len(c["bad_token_count_lines"]["lines"]) < EXAMPLE_CAP:
                c["bad_token_count_lines"]["lines"].append(lineno)
            close_block()  # an unreadable line ends the column's consecutive run
            noncanonical = True
            continue
        c["data_lines"] += 1
        name = t[0]
        if name != cur_col:
            close_block()
            cur_col = name
            if name in col_blocks:
                col_blocks[name] += 1
                c["nonconsecutive_column_blocks"]["count"] += 1
                if len(c["nonconsecutive_column_blocks"]["columns"]) < EXAMPLE_CAP:
                    c["nonconsecutive_column_blocks"]["columns"].append(name)
                noncanonical = True
            else:
                col_blocks[name] = 1
                col_order.append(name)
                if in_int:
                    int_columns.add(name)
                    region_columns.append(name)
        for k in range(1, len(t) - 1, 2):
            row, val = t[k], t[k + 1]
            c["entries_total"] += 1
            kind = row_kind.get(row)
            if kind is None:
                c["unknown_row_references"] += 1
                noncanonical = True
                continue
            if kind == "N":
                c["objective_row_entries"] += 1
                continue
            if kind == "n":
                c["free_n_row_entries_dropped"] += 1
                continue
            c["matrix_entries"] += 1
            if row in cur_block_rows:
                c["duplicate_pair_entries"]["count"] += 1
                if len(c["duplicate_pair_entries"]["examples"]) < EXAMPLE_CAP:
                    c["duplicate_pair_entries"]["examples"].append(
                        {"column": name, "row": row, "line": lineno})
            cur_block_rows.add(row)
    if section == "COLUMNS" and in_int:
        rep["columns"]["markers"]["open_at_columns_section_end"] = True
        open_ended_int_columns = set(region_columns)
    rep["columns"]["distinct"] = len(col_order)
    rep["columns"]["markers"]["columns_marked_integer"] = len(int_columns)
    rep["columns"]["markers"]["columns_marked_integer_open_ended"] = len(open_ended_int_columns)
    if noncanonical:
        rep["tallies_status"] = ("hypothetical: noncanonical constructs present; "
                                 "downstream tallies are for inspection, not reader semantics")
    return rep


def main(argv):
    reports = [scan(p) for p in argv[1:]]
    json.dump(reports[0] if len(reports) == 1 else reports, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
