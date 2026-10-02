#!/usr/bin/env python3
"""Replay tests/repro/*.mps: every reproducer has an expected status in tests/repro/expected.json.

expected.json: { "file.mps": {"expected_status": "optimal|infeasible|unbounded|parse_error|...",
                              "expected_objective": number|null, "known_failing": bool, "category": "...",
                              "source": "case id or hand-written", "observed_bad": "what the engine does",
                              "report": "one paragraph"} }
Outcome per file: OK (engine gives the expected result), KNOWN-FAIL (known_failing and still wrong),
FIXED (known_failing but now correct: remove the flag), REGRESSION (not known_failing and wrong).
Exit code 1 on REGRESSION. The expected status comes from the reference solver (presolve off) and/or from the
construction or the MPS format, never from the engine under test.
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpsio import parse_mps, exact_check, read_sol, OracleError


def run(engine, path, tl):
    js, sol = path + ".json.tmp", path + ".sol.tmp"
    for p in (js, sol):
        if os.path.exists(p):
            os.remove(p)
    try:
        pr = subprocess.run([engine, path, "--time-limit", str(tl), "--json", js, "--sol", sol], capture_output=True, timeout=tl + 30)
    except subprocess.TimeoutExpired:
        return dict(status="harness_timeout", objective=None), None
    try:
        j = json.load(open(js))
    except Exception:
        out = pr.stdout.decode(errors="replace").split()
        st = out[1] if len(out) > 1 and out[0] == "status" else "no_status_rc%s" % pr.returncode
        try:
            ob = float(out[3]) if len(out) > 3 else None
        except ValueError:
            ob = None
        j = dict(status=st, objective=ob, invalid_json=True)
    pt = sol if os.path.exists(sol) else None
    return j, pt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--dir", default="tests/repro")
    ap.add_argument("--time-limit", type=float, default=20.0)
    a = ap.parse_args()
    exp = json.load(open(os.path.join(a.dir, "expected.json"))) if os.path.exists(os.path.join(a.dir, "expected.json")) else {}
    files = sorted(f for f in os.listdir(a.dir) if f.endswith(".mps"))
    rows, regress = [], 0
    for f in files:
        e = exp.get(f)
        if e is None:
            print("repro %s has no entry in expected.json" % f)
            regress += 1
            continue
        p = os.path.join(a.dir, f)
        j, sol = run(os.path.abspath(a.engine), p, a.time_limit)
        ok = j.get("status") == e["expected_status"]
        why = ""
        if ok and e["expected_status"] == "optimal":
            ob = j.get("objective")
            eo = e.get("expected_objective")
            if ob is None or ob != ob or (eo is not None and abs(ob - eo) > 1e-6 * max(1.0, abs(eo))):
                ok, why = False, "objective %r vs expected %r" % (ob, eo)
            elif sol:
                try:
                    chk = exact_check(parse_mps(p), read_sol(sol, parse_mps(p)))
                    if chk["viol_row_rel"] > 1e-6 or chk["viol_bnd_rel"] > 1e-6 or chk["int_viol"] > 1e-6:
                        ok, why = False, "returned point violates the model (row %.2g bound %.2g)" % (chk["viol_row_rel"], chk["viol_bnd_rel"])
                except OracleError as ex:
                    ok, why = False, "oracle cannot read the file: %s" % ex
        elif not ok:
            why = "status %s vs expected %s" % (j.get("status"), e["expected_status"])
        known = e.get("known_failing", False)
        res = "OK" if ok and not known else "FIXED" if ok and known else "KNOWN-FAIL" if known else "REGRESSION"
        if res == "REGRESSION":
            regress += 1
        rows.append((f, e.get("category", ""), e["expected_status"], j.get("status"), res, why))
        for q in (p + ".json.tmp", p + ".sol.tmp"):
            if os.path.exists(q):
                os.remove(q)
    print("%-34s %-16s %-11s %-18s %s" % ("reproducer", "category", "expected", "engine", "outcome"))
    for r in rows:
        print("%-34s %-16s %-11s %-18s %s %s" % (r[0], r[1], r[2], r[3], r[4], ("(" + r[5] + ")") if r[5] else ""))
    from collections import Counter
    c = Counter(r[4] for r in rows)
    print("reproducers: %d  %s" % (len(rows), dict(c)))
    sys.exit(1 if regress else 0)


if __name__ == "__main__":
    main()
