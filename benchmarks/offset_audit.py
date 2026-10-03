#!/usr/bin/env python3
"""Objective-offset audit. Generates tiny LP/MILP/QP files with an objective-row RHS (MPS offset = -RHS), min and max sense,
offsets none/+5.5/-7.25/1e7, and checks the engine objective against HiGHS (reference only) on primal, dual and ipm routes (QP: default and ipm).
Usage: python3 offset_audit.py ./taral  -> prints 'checked N mismatch M'; exit 1 on any mismatch. Needs highspy."""
import glob, json, os, subprocess, sys, tempfile
import highspy
eng = sys.argv[1]; d = tempfile.mkdtemp()
def lp(name, sense, off, integer=False, quad=""):
    a, b = (" M0 'MARKER' 'INTORG'\n", " M1 'MARKER' 'INTEND'\n") if integer else ("", "")
    s = f"NAME {name}\n" + ("OBJSENSE\n MAX\n" if sense == "max" else "") + "ROWS\n N OBJ\n G R1\n L R2\nCOLUMNS\n" + a
    s += " X OBJ 1 R1 1\n X R2 1\n Y OBJ 2 R1 1\n Y R2 2\n" + b + "RHS\n RHS R1 3\n RHS R2 10\n" + (f" RHS OBJ {off}\n" if off is not None else "")
    open(os.path.join(d, name + ".mps"), "w").write(s + "BOUNDS\n UP B X 8\n UP B Y 8\n" + quad + "ENDATA\n")
for sn in ("min", "max"):
    for off, tag in ((None, "none"), (5.5, "pos"), (-7.25, "neg"), (1e7, "big")):
        lp(f"o_{sn}_{tag}", sn, off); lp(f"om_{sn}_{tag}", sn, off, True)
        lp(f"oq_{sn}_{tag}", sn, off, quad="QUADOBJ\n X X " + ("-" if sn == "max" else "") + "2\n")
bad = n = 0
for f in sorted(glob.glob(d + "/o*.mps")):
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.readModel(f); h.run(); ref = h.getInfo().objective_function_value
    for m in ([None, "ipm"] if "/oq_" in f else [None, "dual", "ipm"]):
        js = os.path.join(d, "j.json")
        subprocess.run([eng, f, "--json", js] + (["--method", m] if m else []), capture_output=True)
        r = json.load(open(js)); n += 1
        if not (r.get("status") == "optimal" and abs(r["objective"] - ref) <= 1e-6 * (1 + abs(ref))):
            bad += 1; print("MISMATCH", os.path.basename(f), m, r.get("status"), r.get("objective"), ref)
print("checked", n, "mismatch", bad); sys.exit(1 if bad else 0)
