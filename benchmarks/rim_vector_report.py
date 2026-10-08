#!/usr/bin/env python3
"""Static MPS rim-vector report. Usage: rim_vector_report.py file.mps [file2.mps ...]

Lists, for the RHS / RANGES / BOUNDS sections of an MPS file, which named rim-vector
sets are present, which name a first-set reader selects (the first non-empty set name
encountered, the documented CPLEX/IMSL default policy), and how many entries a
first-set reader would skip because they carry a different set name. Nameless entries
(no set name) are counted separately: under TARAL's free-token policy, nameless
entries are counted as applied. That nameless rule is TARAL-specific, not documented
for all first-set readers.

This is a STATIC checker: it does not run or modify any solver and it assumes the
free-format tokenization below. It is not a general parser, and it diverges from the
reader on invalid or noncanonical input. A repeated RHS / RANGES / BOUNDS section
header is noncanonical; occurrences are listed separately and NO selection semantics
are assumed for that kind ("semantics": "not_assumed"). The per-occurrence
selection/skip tallies are local hypotheticals for inspection, not model semantics.
Scanning stops at ENDATA; any trailing content is ignored and not compared against
the reader. Unknown column-1 headers and malformed data lines are counted, not
interpreted.
"""
import hashlib, json, sys

SECTIONS = {"NAME", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "OBJSENSE",
            "OBJNAME", "QUADOBJ", "QMATRIX", "QSECTION", "SOS", "ENDATA"}
RIM = ("RHS", "RANGES", "BOUNDS")
NO_VALUE_BOUNDS = {"FR", "MI", "PL", "BV"}


def bound_set_name(t):
    """Set name of a BOUNDS data line under free-format rules (mirrors the reader):
    4 tokens -> t[1]; 3 tokens -> t[1] for no-value types, else none; 2 tokens -> none."""
    if len(t) == 4:
        return t[1]
    if len(t) > 4:
        return None  # malformed: too many fields
    if len(t) == 3:
        return t[1] if t[0] in NO_VALUE_BOUNDS else ""
    if len(t) == 2:
        return ""
    return None  # malformed


def scan(path):
    with open(path, "rb") as f:
        raw = f.read()
    report = {"file": path, "sha256": hashlib.sha256(raw).hexdigest(),
              "policy": "first-set (first non-empty set name wins; nameless entries always applied)",
              "rhs": None, "ranges": None, "bounds": None}
    kind_state = {k: {"occurrences": [], "repeated": False} for k in RIM}
    unknown_headers = []
    section, lineno = None, 0
    for line in raw.decode("utf-8", "replace").splitlines():
        lineno += 1
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        t = s.split()
        head = t[0].upper()
        # Section headers start in column 1 (MPS convention); indented data lines
        # whose first token happens to equal a keyword are not headers.
        if line[0] not in " \t" and head in SECTIONS:
            if head == "ENDATA":
                break
            section = head
            if section in RIM:
                occ = {"occurrence": len(kind_state[section]["occurrences"]) + 1,
                       "header_line": lineno, "names": [], "selected": None,
                       "entries_by_name": {}, "name_skipped": 0, "nameless": 0}
                if len(kind_state[section]["occurrences"]) >= 1:
                    kind_state[section]["repeated"] = True
                kind_state[section]["occurrences"].append(occ)
            continue
        if line[0] not in " \t":
            section = None
            unknown_headers.append(lineno)
            continue
        if section not in RIM:
            continue
        occ = kind_state[section]["occurrences"][-1]
        if section in ("RHS", "RANGES"):
            start = 1 if len(t) % 2 == 1 else 0  # odd token count: leading set name
            name = t[0] if start == 1 else ""
            pairs = (len(t) - start) // 2
            if pairs == 0:
                occ.setdefault("malformed_lines", 0)
                occ["malformed_lines"] += 1
                continue
        else:  # BOUNDS
            name = bound_set_name(t)
            if name is None:
                occ.setdefault("malformed_lines", 0)
                occ["malformed_lines"] += 1
                continue
            pairs = 1
        if name == "":
            occ["nameless"] += pairs
            continue
        if not occ["names"]:
            occ["selected"] = name
        if name not in occ["names"]:
            occ["names"].append(name)
        occ["entries_by_name"][name] = occ["entries_by_name"].get(name, 0) + pairs
        if name != occ["selected"]:
            occ["name_skipped"] += pairs
    if unknown_headers:
        report["unknown_column1_headers_at_lines"] = unknown_headers
    for kind in RIM:
        st = kind_state[kind]
        if not st["occurrences"]:
            continue
        out = {"occurrences": st["occurrences"], "repeated": st["repeated"]}
        if st["repeated"]:
            out["semantics"] = "not_assumed"
            out["note"] = ("repeated " + kind + " section header (noncanonical); "
                           "occurrences listed separately, no selection semantics assumed")
        report[kind.lower()] = out
    return report


def main():
    reports = [scan(p) for p in sys.argv[1:]]
    print(json.dumps(reports if len(reports) != 1 else reports[0], indent=2))


if __name__ == "__main__":
    main()
