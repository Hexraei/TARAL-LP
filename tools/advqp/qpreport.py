#!/usr/bin/env python3
"""Per-category QP table from results.jsonl, with absolute violations and KKT residuals. Exit 1 on any NEW failure.

A failing case is KNOWN when its id is listed in tests/qp_known_failures.tsv (id <TAB> class <TAB> reproducer <TAB> detail).
Columns: pass, pass* (engine certified, reference/construction inexact), partial (dual_infeasible), border (Q within
noise of semidefinite), wrong, noans, other. pass-rate = (pass+pass*) / cases; strict = same but also requiring the
exact ABSOLUTE row/bound violation <= 1e-6 (the gate itself is the 1+|bound| relative one).
"""
import argparse
import json
import os
import sys
from collections import Counter, OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qpcase

WRONG = ("FAIL_STATUS", "FAIL_OBJ", "FAIL_FEAS", "FAIL_KKT", "FAIL_NOSOL", "PARSE_ERR")
OTHER = ("REF_UNRESOLVED", "GEN_MISMATCH", "HARNESS_ERR")


def load_known(path):
    known = {}
    if os.path.exists(path):
        for ln in open(path):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            known[p[0]] = p[1] if len(p) > 1 else ""
    return known


def absv(r):
    a = [r.get("row_abs"), r.get("bnd_abs")]
    a = [v for v in a if v is not None]
    return max(a) if a else None


def pct(a, b):
    return "%.1f%%" % (100.0 * a / b) if b else "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--known", default="tests/qp_known_failures.tsv")
    ap.add_argument("--md", default="")
    a = ap.parse_args()
    recs = [json.loads(l) for l in open(a.results)]
    known = load_known(a.known)
    cats = OrderedDict((c, Counter()) for c, _, _ in qpcase.ALL if any(r["cat"] == c for r in recs))
    mx = {c: dict(abs=0.0, stat=0.0, gap=0.0, act=0.0, xerr=0.0, n_opt=0) for c in cats}
    new_fail = []
    for r in recs:
        v, c = r["verdict"], cats[r["cat"]]
        c["n"] += 1
        key = {"PASS": "pass", "PASS_REF": "ref", "PARTIAL": "partial", "BORDER": "border", "NOANS": "noans"}.get(v)
        if key is None:
            key = "wrong" if v in WRONG else "other"
        c[key] += 1
        av = absv(r)
        if av is not None:
            mx[r["cat"]]["n_opt"] += 1
            mx[r["cat"]]["abs"] = max(mx[r["cat"]]["abs"], av)
            if av > 1e-6:
                c["abs_gt"] += 1
            if max(r.get("row_rel") or 0, r.get("bnd_rel") or 0) > 1e-6:
                c["rel_gt"] += 1
            if r.get("stat_rel") is not None:
                mx[r["cat"]]["stat"] = max(mx[r["cat"]]["stat"], r["stat_rel"])
                mx[r["cat"]]["gap"] = max(mx[r["cat"]]["gap"], r.get("gap_rel") or 0)
                mx[r["cat"]]["act"] = max(mx[r["cat"]]["act"], r.get("act_rel") or 0)
                mx[r["cat"]]["xerr"] = max(mx[r["cat"]]["xerr"], r.get("x_err") or 0)
                if (r.get("act_rel") or 0) > 1e-6:
                    c["act_gt"] += 1
            if v in ("PASS", "PASS_REF") and av > 1e-6:
                c["strict_loss"] += 1
        if key in ("wrong", "noans", "other"):
            if r["id"] in known:
                c["known"] += 1
            else:
                new_fail.append(r)
    hdr = "%-22s %5s %5s %5s %7s %6s %5s %5s %5s %6s %7s %7s" % (
        "category", "cases", "pass", "pass*", "partial", "border", "wrong", "noans", "other", "known", "rate", "strict")
    lines = [hdr, "-" * len(hdr)]
    tot = Counter()
    for cat, c in cats.items():
        ok = c["pass"] + c["ref"]
        strict = ok - c["strict_loss"]
        lines.append("%-22s %5d %5d %5d %7d %6d %5d %5d %5d %6d %7s %7s" % (
            cat, c["n"], c["pass"], c["ref"], c["partial"], c["border"], c["wrong"], c["noans"], c["other"], c["known"],
            pct(ok, c["n"]), pct(strict, c["n"])))
        tot.update(c)
    ok = tot["pass"] + tot["ref"]
    lines.append("-" * len(hdr))
    lines.append("%-22s %5d %5d %5d %7d %6d %5d %5d %5d %6d %7s %7s" % (
        "TOTAL", tot["n"], tot["pass"], tot["ref"], tot["partial"], tot["border"], tot["wrong"], tot["noans"], tot["other"],
        tot["known"], pct(ok, tot["n"]), pct(ok - tot["strict_loss"], tot["n"])))
    table = "\n".join(lines)
    print(table)
    # violation / certificate table (only cases where the engine returned a point)
    l2 = ["%-22s %6s %12s %9s %9s %12s %12s %12s %9s %10s" % ("category", "points", "max abs viol", "abs>1e-6", "rel>1e-6",
                                                              "max stat rel", "max gap rel", "max act rel", "act>1e-6", "max |x-x*|")]
    l2.append("-" * len(l2[0]))
    for cat, c in cats.items():
        m = mx[cat]
        l2.append("%-22s %6d %12.3e %9d %9d %12.3e %12.3e %12.3e %9d %10.2e" % (
            cat, m["n_opt"], m["abs"], c["abs_gt"], c["rel_gt"], m["stat"], m["gap"], m["act"], c["act_gt"], m["xerr"]))
    print()
    print("\n".join(l2))
    # failures by class
    cls = Counter()
    for r in recs:
        v = r["verdict"]
        if v in WRONG or v in OTHER or v == "NOANS":
            cls[(r["cat"], v, r.get("ours_status"))] += 1
    print()
    print("%-22s %-14s %-18s %s" % ("category", "verdict", "engine status", "cases"))
    for (cat, v, st), n in sorted(cls.items()):
        print("%-22s %-14s %-18s %d" % (cat, v, st, n))
    print()
    print("known: %d   NEW failures: %d" % (tot["known"], len(new_fail)))
    for r in new_fail[:40]:
        print("  NEW %s %s: %s" % (r["id"], r["verdict"], (r.get("detail") or "")[:150]))
    if a.md:
        with open(a.md, "w") as f:
            f.write("```\n" + table + "\n\n" + "\n".join(l2) + "\n```\n")
    sys.exit(1 if new_fail else 0)


if __name__ == "__main__":
    main()
