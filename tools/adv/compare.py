#!/usr/bin/env python3
"""Per-category before/after table across solver paths, from harness results.jsonl files.

  compare.py --before simplex=out/methods/simplex/results.jsonl --after dual=.../results.jsonl ipm=.../results.jsonl
             [--known tests/known_failures.tsv] [--md out.md]

'before' is the default primal-simplex path, 'after' the other paths, run on the same cases with the same gates. Per
category and path: cases, ok (pass + pass*), wrong (wrong status/objective/infeasible point), noans (time/iteration limit,
numerical_failure, dual_infeasible and similar non-answers), other (reference unresolved, harness error).
A failure of a path is
  shared   - the same case also fails on the before path (a known simplex failure when listed in known_failures.tsv)
  NEW      - the case passes on the before path and fails on this one (a path-specific finding)
and a case that fails before but passes on this path is reported as fixed. Exit status 1 when any path has a NEW failure
(or a failing case that the before path never ran), 0 otherwise; the NEW cases are listed with their detail.
"""
import argparse
import json
import os
import sys
from collections import Counter, OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from report import WRONG, load_known

ORDER = [c for c, _, _ in gen.ALL + gen.BIG_CATS]


def load(path):
    return {r["id"]: r for r in (json.loads(l) for l in open(path))}


def kind(v):
    if v in ("PASS", "PASS_REF"):
        return "ok"
    return "wrong" if v in WRONG else "noans" if v == "NOANS" else "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True, help="NAME=results.jsonl")
    ap.add_argument("--after", nargs="+", required=True, help="NAME=results.jsonl ...")
    ap.add_argument("--known", default="tests/known_failures.tsv")
    ap.add_argument("--md", default="")
    a = ap.parse_args()
    bname, bpath = a.before.split("=", 1)
    paths = OrderedDict([(bname, load(bpath))])
    for s in a.after:
        n, p = s.split("=", 1)
        paths[n] = load(p)
    names = list(paths)
    known = load_known(a.known)
    base = paths[bname]
    cats = [c for c in ORDER if any(r["cat"] == c for r in base.values())]
    ids_by_cat = {c: sorted(i for i, r in base.items() if r["cat"] == c) for c in cats}

    def cell(name, c):
        cnt = Counter()
        for i in ids_by_cat[c]:
            r = paths[name].get(i)
            if r is None:
                cnt["not run"] += 1
            else:
                cnt[kind(r["verdict"])] += 1
        return cnt

    lines = []
    head = "%-17s %-5s %5s |" % ("category", "kind", "cases")
    for n in names:
        head += " %-26s|" % (n + " ok/wrong/noans/other")
    lines += [head, "-" * len(head)]
    tot = {n: {k: Counter() for k in ("LP", "MILP")} for n in names}
    for c in cats:
        k = "MILP" if gen.IS_MIP.get(c) or any(r.get("mip") for r in base.values() if r["cat"] == c) else "LP"
        row = "%-17s %-5s %5d |" % (c, k, len(ids_by_cat[c]))
        for n in names:
            cnt = cell(n, c)
            tot[n][k].update(cnt)
            if cnt["not run"]:
                row += " %-26s|" % "not run (method n/a)"
            else:
                n_ = sum(cnt.values())
                row += " %4d/%3d/%3d/%3d %6.1f%%  |" % (cnt["ok"], cnt["wrong"], cnt["noans"], cnt["other"], 100.0 * cnt["ok"] / max(1, n_))
        lines.append(row)
    lines.append("-" * len(head))
    for k in ("LP", "MILP"):
        n_cases = sum(sum(tot[bname][k].values()) for _ in [0])
        if not n_cases:
            continue
        row = "%-17s %-5s %5d |" % ("TOTAL", k, n_cases)
        for n in names:
            cnt = tot[n][k]
            if cnt["not run"]:
                row += " %-26s|" % "not run (method n/a)"
            else:
                row += " %4d/%3d/%3d/%3d %6.1f%%  |" % (cnt["ok"], cnt["wrong"], cnt["noans"], cnt["other"], 100.0 * cnt["ok"] / max(1, sum(cnt.values())))
        lines.append(row)

    findings = []
    nnew = 0
    for n in names[1:]:
        new, shared, fixed, missing = [], [], [], []
        for i, r in sorted(paths[n].items()):
            if r["cat"] not in cats:
                continue
            bk = base.get(i)
            if bk is None:
                missing.append(r)
                continue
            bf, af = kind(bk["verdict"]) != "ok", kind(r["verdict"]) != "ok"
            if af and bf:
                shared.append(r)
            elif af:
                new.append(r)
            elif bf:
                fixed.append(bk)
        nnew += len(new) + len(missing)
        findings.append("")
        findings.append("== %s vs %s: %d failing (%d shared with %s, %d NEW), %d fixed (fail on %s, pass here)" % (
            n, bname, len(shared) + len(new), len(shared), bname, len(new), len(fixed), bname))
        by = Counter((r["cat"], r["verdict"]) for r in new)
        for (cat, v), cnt in sorted(by.items()):
            findings.append("   NEW %-18s %-12s %d" % (cat, v, cnt))
        for r in new[:80]:
            findings.append("   NEW %-24s %-11s %s | %s" % (r["id"], r["verdict"], r.get("ours_status"), (r.get("detail") or "")[:130]))
        if len(new) > 80:
            findings.append("   ... %d more NEW" % (len(new) - 80))
        if fixed:
            findings.append("   fixed here: " + ", ".join(sorted(f["id"] for f in fixed)[:40]) + (" ..." if len(fixed) > 40 else ""))
        for r in missing[:10]:
            findings.append("   NOT RUN on %s: %s" % (bname, r["id"]))
    out = "\n".join(lines + findings)
    print(out)
    if a.md:
        with open(a.md, "w") as f:
            f.write("```\n" + out + "\n```\n")
    sys.exit(1 if nnew else 0)


if __name__ == "__main__":
    main()
