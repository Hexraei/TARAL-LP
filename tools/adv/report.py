#!/usr/bin/env python3
"""Pass/fail table by category from harness results (+ parser table + repro replay). Exit 1 on any NEW failure.

A failing case is KNOWN when its id is listed in tests/known_failures.tsv (id <TAB> class <TAB> note); the file
is the committed record of failures that have a reproducer under tests/repro/.
"""
import argparse
import json
import os
import sys
from collections import Counter, OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen

WRONG = ("FAIL_STATUS", "FAIL_OBJ", "FAIL_FEAS", "FAIL_NOSOL")
OTHER = ("PARSE_ERR", "REF_UNRESOLVED", "GEN_MISMATCH", "ORACLE_ERR", "HARNESS_ERR")


def load_known(path):
    known = {}
    if os.path.exists(path):
        for ln in open(path):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            known[p[0]] = p[1] if len(p) > 1 else ""
    return known


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--known", default="tests/known_failures.tsv")
    ap.add_argument("--md", default="")
    a = ap.parse_args()
    recs = [json.loads(l) for l in open(a.results)]
    known = load_known(a.known)
    cats = OrderedDict((c, Counter()) for c, _, _ in gen.ALL + gen.BIG_CATS if any(r["cat"] == c for r in recs))
    new_fail = []
    for r in recs:
        v = r["verdict"]
        c = cats[r["cat"]]
        c["n"] += 1
        if v == "PASS":
            c["pass"] += 1
        elif v == "PASS_REF":
            c["ref"] += 1
        else:
            kind = "wrong" if v in WRONG else "noans" if v == "NOANS" else "other"
            c[kind] += 1
            c["known" if r["id"] in known else "new"] += 1
            if r["id"] not in known:
                new_fail.append(r)
    head = "%-18s %-5s %6s %6s %7s %6s %6s %6s %6s %8s" % ("category", "kind", "cases", "pass", "pass*", "wrong", "noans", "other", "known", "pass-rate")
    lines = [head, "-" * len(head)]
    tot = {k: Counter() for k in ("LP", "MILP")}
    for cat, c in cats.items():
        kind = "MILP" if gen.IS_MIP.get(cat) else "LP"
        ok = c["pass"] + c["ref"]
        lines.append("%-18s %-5s %6d %6d %7d %6d %6d %6d %6d %7.1f%%" % (cat, kind, c["n"], c["pass"], c["ref"], c["wrong"], c["noans"], c["other"],
                                                                         c["known"], 100.0 * ok / max(1, c["n"])))
        tot[kind].update(c)
    lines.append("-" * len(head))
    for kind in ("LP", "MILP"):
        c = tot[kind]
        ok = c["pass"] + c["ref"]
        lines.append("%-18s %-5s %6d %6d %7d %6d %6d %6d %6d %7.1f%%" % ("TOTAL", kind, c["n"], c["pass"], c["ref"], c["wrong"], c["noans"], c["other"],
                                                                         c["known"], 100.0 * ok / max(1, c["n"])))
    lines.append("pass* = engine matched after the reference was shown wrong/inexact (presolve-off re-solve, exact check, or construction)")
    lines.append("wrong = wrong status/objective/feasibility; noans = time/iteration limit, numerical_failure, unresolved MILP; other = harness/generator")
    out = "\n".join(lines)
    print(out)
    if a.md:
        with open(a.md, "w") as f:
            f.write("```\n" + out + "\n```\n")
    if new_fail:
        print("\nNEW (not in known_failures.tsv): %d" % len(new_fail))
        for r in new_fail[:60]:
            print("  %-26s %-12s %s" % (r["id"], r["verdict"], (r.get("detail") or "")[:150]))
    sys.exit(1 if new_fail else 0)


if __name__ == "__main__":
    main()
