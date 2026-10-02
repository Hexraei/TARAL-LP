#!/usr/bin/env python3
"""Adversarial MILPs whose node LPs carry roundoff stationarity on columns without a finite bound.
Each must solve to optimal with the recorded objective; an unverified node LP must never turn a
solvable model into numerical_failure. Usage: node_lp_proof_regression.py --engine out/taral"""
import argparse, json, os, subprocess, sys, tempfile
EXPECT = {"int_free_neg-51008": (-107.25, 9), "int_free_neg-51009": (93.0, 3),
          "markers_mixed-52000": (-191.911886182, 343), "markers_mixed-52005": (263.081433264, 155),
          "objconst_max-54001": (495.733333333, 7470)}
ap = argparse.ArgumentParser(); ap.add_argument("--engine", required=True); a = ap.parse_args()
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "node_lp_proof")
bad = 0
for name, (obj, nodes) in EXPECT.items():
    with tempfile.TemporaryDirectory() as t:
        js = os.path.join(t, "r.json")
        subprocess.run([a.engine, os.path.join(d, name + ".mps"), "--time-limit", "60", "--json", js],
                       check=False, capture_output=True)
        r = json.load(open(js))
    ok = (r.get("status") == "optimal" and abs(r["objective"] - obj) <= 1e-6 * max(1, abs(obj))
          and abs(r["best_bound"] - obj) <= 1e-6 * max(1, abs(obj)))
    print(("PASS " if ok else "FAIL ") + name, r.get("status"), r.get("objective"), r.get("best_bound"), r.get("nodes"))
    bad += not ok
sys.exit(1 if bad else 0)
