"""Regression: a badly scaled LP on which the plain solve ends in numerical_failure and the KKT-gated
equilibrated retry recovers the HiGHS optimum. Usage: python tests/scaled_retry_check.py --engine ./taral"""
import argparse, subprocess, sys, os
ap = argparse.ArgumentParser(); ap.add_argument("--engine", required=True); a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([a.engine, os.path.join(here, "fixtures", "scaled_retry_c87.mps"), "--time-limit", "20"],
                     capture_output=True, text=True).stdout.strip().split("\n")[-1]
t = out.split()
ok = t[1] == "optimal" and abs(float(t[3]) - (-7350579.038607513)) <= 1e-6 * 7350579.04
print(("PASS " if ok else "FAIL ") + out[:120]); sys.exit(0 if ok else 1)
