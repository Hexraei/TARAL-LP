#!/usr/bin/env python3
"""Second reference for the Maros-Meszaros rows that have no usable default-HiGHS reference.

For every ledger row whose HiGHS reference is not Optimal, solve the same file with HiGHS and
qp_regularization_value=1e-12 (the setting that made HiGHS agree with taral on QBORE3D, GOULDQP2, PRIMALC8) and compare
with taral's ledger objective under the ledger rule (|diff| <= 1e-6 * max(1, |ref|)). taral's row/bound violations come
from the ledger's independent recheck. A row stays "unreferenced" if HiGHS again returns no optimal answer.

Usage: ~/.venvs/taral-gpu/bin/python benchmarks/qp_recheck_unreferenced.py [--highs-limit 300] [--workers 4]
Writes results/qp_maros_meszaros/recheck_unreferenced.csv.
"""
import argparse, csv, os, sys, time
from multiprocessing import Pool

import highspy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qp_maros_meszaros as mm

LEDGER = mm.OUT
OUT = os.path.join(os.path.dirname(LEDGER), "recheck_unreferenced.csv")


def solve(args):
    name, path, limit = args
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(limit))
    h.setOptionValue("qp_regularization_value", 1e-12)
    h.readModel(path)
    t = time.time()
    try:
        h.run()
    except Exception as e:
        return name, "error:" + str(e)[:60], None, time.time() - t
    ok = h.getModelStatus() == highspy.HighsModelStatus.kOptimal
    return name, h.modelStatusToString(h.getModelStatus()), h.getInfo().objective_function_value if ok else None, time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--highs-limit", type=float, default=300)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    ledger = {r["instance"]: r for r in csv.DictReader(open(LEDGER))}
    todo = [n for n, r in ledger.items() if r["ref_status"] != "Optimal"]
    src = {n: (s, u) for n, s, u in mm.instances(10 ** 9)}
    jobs = [(n, mm.fetch(n, *src[n]), a.highs_limit) for n in todo]
    with Pool(a.workers) as p:
        res = p.map(solve, jobs, chunksize=1)
    rows = []
    for name, st, ref, wall in sorted(res):
        r = ledger[name]
        obj = float(r["objective"]) if r["objective"] else None
        viol = max(float(r["chk_row_viol"] or 0), float(r["chk_bound_viol"] or 0))
        err = mm.rel_err(obj, ref)
        if ref is None:
            verdict = "unreferenced"
        elif r["status"] == "optimal" and err <= mm.OBJ_TOL and viol <= mm.VIOL_TOL:
            verdict = "agree"
        else:
            verdict = "DISAGREE"
        rows.append([name, r["rows"], r["cols"], r["ref_status"], st, r["status"], obj, ref, err, viol, round(wall, 1), verdict])
        print(*rows[-1], flush=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance", "rows", "cols", "default_highs", "highs_reg1e-12", "taral_status", "taral_objective",
                    "highs_reg_objective", "rel_err", "taral_max_violation", "highs_wall_s", "verdict"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
