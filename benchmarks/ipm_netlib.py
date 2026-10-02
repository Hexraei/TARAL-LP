"""IPM line on the pinned Netlib corpus (separate from the simplex headline).
Usage: python benchmarks/ipm_netlib.py [--workers 4] [--limit 60] [case ...]
Writes out/ipm/netlib_ipm.csv. Reference objectives: HiGHS on the original files (results/cpp_818ffec_60s/highs_reference.json).
Objective rel error = |ipm - highs| / max(1, |highs|), objective recomputed independently by recheck().
"""
import argparse, csv, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ipm_check import ROOT, run_ipm, recheck, rel_err

ap = argparse.ArgumentParser()
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--limit", type=float, default=60.0)
ap.add_argument("cases", nargs="*")
args = ap.parse_args()
ref = {k: {"status": v["status"], "objective": v["objective"]}
       for k, v in json.loads((ROOT / "results" / "cpp_818ffec_60s" / "highs_reference.json").read_text()).items()}
files = sorted((ROOT / "corpus").glob("*.mps"))
if args.cases:
    files = [f for f in files if f.stem in args.cases]
work = ROOT / "out" / "ipm" / "runs"


def one(f):
    r, sol = run_ipm(f, work, args.limit)
    row = {"case": f.stem, "status": r.get("status"), "iterations": r.get("iterations"), "wall_s": r.get("wall_s"),
           "objective": r.get("objective"), "primal_res": r.get("primal_res"), "dual_res": r.get("dual_res"),
           "gap": r.get("gap"), "max_row_viol": r.get("max_row_viol"), "kkt_dim": r.get("kkt_dim"),
           "factor_nnz": r.get("factor_nnz"), "message": r.get("message", "")}
    h = ref.get(f.stem, {})
    row["highs_status"] = h.get("status")
    row["highs_obj"] = h.get("objective")
    if sol:
        row.update(recheck(f, sol))
    row["obj_rel_err"] = rel_err(row.get("chk_obj"), row["highs_obj"])
    return row


with ThreadPoolExecutor(args.workers) as ex:
    rows = list(ex.map(one, files))
out = ROOT / "out" / "ipm" / "netlib_ipm.csv"
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, keys)
    w.writeheader()
    w.writerows(rows)
for r in rows:
    e = r["obj_rel_err"]
    print(f"{r['case']:12s} {r['status']:18s} it {r['iterations']!s:>4} {r['wall_s'] or 0:8.2f}s "
          f"relerr {e if e is None else f'{e:.1e}'} pres {r['primal_res']} dres {r['dual_res']} gap {r['gap']} {r['message']}")
ok = [r for r in rows if r["status"] == "optimal" and r["case"] != "truss"]
print(f"optimal (93 Netlib, truss excluded): {len(ok)}/{sum(r['case'] != 'truss' for r in rows)}; wrote {out}")
