#!/usr/bin/env python3
"""Install a reproducer into tests/repro/ with its expected result and report paragraph.

Expected status/objective come from the reference, never from the engine: HiGHS 1.15.1 is run three ways
(simplex presolve off, simplex presolve on, interior point presolve off). They must agree (status, and objective to
1e-6 relative) and the returned point must satisfy the original model in exact arithmetic (<= 1e-9 relative). If
they do not, the install is refused unless --expected-status is given (hand-derived truth, e.g. from the MPS format
or an exact construction), in which case --why-expected records the derivation.

Usage: mkrepro.py --src FILE.mps --name NAME --category CAT --source "cat-seed|hand-written" --observed "..." \
                  --report "one paragraph" [--known/--not-known] [--engine ./taral]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpsio import parse_mps, exact_check
from ref import solve_ref


def reference_truth(path):
    m = parse_mps(path)
    runs = [("simplex/presolve-off", solve_ref(m, presolve=False, time_limit=30)),
            ("simplex/presolve-on", solve_ref(m, presolve=True, time_limit=30)),
            ("ipm/presolve-off", solve_ref(m, presolve=False, time_limit=30, solver="ipm"))]
    det = [(n, r) for n, r in runs if r["status"] in ("optimal", "infeasible", "unbounded")]
    sts = {r["status"] for _, r in det}
    if len(sts) != 1 or len(det) < 2:
        return None, "reference configurations disagree: %s" % [(n, r["status"]) for n, r in runs]
    st = sts.pop()
    obj = None
    if st == "optimal":
        objs = [r["objective"] for _, r in det]
        if max(objs) - min(objs) > 1e-6 * max(1.0, abs(objs[0])):
            return None, "reference objectives disagree: %s" % objs
        obj = objs[0]
        for n, r in det:
            c = exact_check(m, r["x"])
            if c["viol_row_rel"] > 1e-9 or c["viol_bnd_rel"] > 1e-9 or c["int_viol"] > 1e-9:
                return None, "reference point from %s is not exactly feasible (row %.2g bound %.2g)" % (n, c["viol_row_rel"], c["viol_bnd_rel"])
    return (st, obj), "HiGHS 1.15.1: %s agree" % ", ".join(n for n, _ in det)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--category", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--observed", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--known", dest="known", action="store_true", default=True)
    ap.add_argument("--not-known", dest="known", action="store_false")
    ap.add_argument("--expected-status")
    ap.add_argument("--expected-objective", type=float)
    ap.add_argument("--why-expected", default="")
    ap.add_argument("--dir", default="tests/repro")
    ap.add_argument("--method", default="", help="engine --method the failure belongs to (dual|ipm); empty = default path")
    ap.add_argument("--time-limit", type=float, default=0, help="per-reproducer engine time limit for the replay (0 = replay default)")
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    dst = os.path.join(a.dir, a.name + ".mps")
    if a.expected_status:
        st, obj, why = a.expected_status, a.expected_objective, a.why_expected
    else:
        t, why = reference_truth(a.src)
        if t is None:
            print("REFUSED:", why)
            sys.exit(1)
        st, obj = t
    shutil.copyfile(a.src, dst)
    path = os.path.join(a.dir, "expected.json")
    exp = json.load(open(path)) if os.path.exists(path) else {}
    exp[a.name + ".mps"] = dict(category=a.category, source=a.source, expected_status=st, expected_objective=obj,
                                expected_basis=why, known_failing=a.known, observed_bad=a.observed, report=a.report)
    if a.method:
        exp[a.name + ".mps"]["method"] = a.method
    if a.time_limit:
        exp[a.name + ".mps"]["time_limit"] = a.time_limit
    json.dump(dict(sorted(exp.items())), open(path, "w"), indent=1)
    print("installed %s: expected %s %s (%s)" % (dst, st, obj, why))


if __name__ == "__main__":
    main()
