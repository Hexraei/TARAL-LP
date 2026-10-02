#!/usr/bin/env python3
"""Replay tests/qp_repro/*.mps. expected.json maps each file to {"case": "<category>-<seed>" or null, ...}.

A reproducer that names a generator case is re-judged with the full harness (exact violations, KKT certificate, HiGHS,
construction) after checking that the file is byte-identical to what the generator writes today; a hand-written
reproducer (case null) is judged on status and objective only. Outcome per file:
  OK          the engine now gives an acceptable result (verdict PASS / PASS_REF / PARTIAL when partial_ok)
  KNOWN-FAIL  known_failing and still wrong (the verdict is printed)
  FIXED       known_failing but now acceptable: remove the flag
  REGRESSION  not known_failing and wrong, or the file no longer matches its generator case
Exit code 1 on REGRESSION.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qpcase
import qpharness

GOOD = ("PASS", "PASS_REF", "BORDER")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--dir", default="tests/qp_repro")
    ap.add_argument("--time-limit", type=float, default=30.0)
    a = ap.parse_args()
    exp = json.load(open(os.path.join(a.dir, "expected.json")))
    regress = 0
    for f in sorted(x for x in os.listdir(a.dir) if x.endswith(".mps")):
        e = exp.get(f)
        if e is None:
            print("%-44s no entry in expected.json" % f)
            regress += 1
            continue
        path = os.path.join(a.dir, f)
        sol, js = path + ".sol.tmp", path + ".json.tmp"
        ours = qpharness.run_engine(os.path.abspath(a.engine), path, sol, js, a.time_limit)
        verdict, detail = None, ""
        if e.get("case"):
            cat, seed = e["case"].rsplit("-", 1)
            case = qpcase.gen(cat, int(seed) - qpcase.BASE[cat])
            if qpcase.write_mps(case) != open(path).read():
                print("%-44s REGRESSION  file differs from generator case %s" % (f, e["case"]))
                regress += 1
                continue
            x = qpharness.read_sol(sol, case.n) if ours.get("status") == "optimal" else None
            rec = qpharness.judge(case, ours, x, a.time_limit)
            verdict, detail = rec["verdict"], rec.get("detail", "")
        else:
            ok = ours.get("status") == e["expected_status"]
            if ok and e.get("expected_objective") is not None and ours.get("objective") is not None:
                ok = abs(ours["objective"] - e["expected_objective"]) <= 1e-6 * max(1.0, abs(e["expected_objective"]))
            verdict = "PASS" if ok else "FAIL_STATUS"
            detail = "engine %s %s" % (ours.get("status"), ours.get("objective"))
        for p in (sol, js):
            if os.path.exists(p):
                os.remove(p)
        good = verdict in GOOD or (verdict == "PARTIAL" and e.get("partial_ok"))
        if good and e.get("known_failing"):
            out = "FIXED"
        elif good:
            out = "OK"
        elif e.get("known_failing"):
            out = "KNOWN-FAIL"
        else:
            out = "REGRESSION"
            regress += 1
        print("%-44s %-10s %s %s" % (f, out, verdict, detail[:110]))
    sys.exit(1 if regress else 0)


if __name__ == "__main__":
    main()
