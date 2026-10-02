#!/usr/bin/env python3
"""Maros-Meszaros convex-QP ledger for the C++ interior-point QP solver (local measurement, not a published result).

Source: the original QPS files of the Maros-Meszaros set (138 problems) as mirrored in
github.com/YimingYAN/QP-Test-Problems, folder QPS_Files, pinned to commit 871cc300 and by the sha256 of every
file (results/qp_maros_meszaros/sources.csv; files are downloaded on demand, never committed). Only instances with
rows + cols <= --max-size (default 5000) are run. Reference objective: a HiGHS solve of the SAME file (cached in
out/maros_meszaros/highs_ref.json). Every taral solution is re-checked against the ORIGINAL model as read by
HiGHS (rows, bounds, objective c'x + 0.5 x'Qx + offset), so the check shares no code with taral's reader.

Convention: QUADOBJ lists one triangle of Q and the objective is 0.5 x'Qx (checked on a hand case first:
min x + 0.5*2*x^2, x >= 1 gives 2; an off-diagonal entry listed once is mirrored). All 138 files use QUADOBJ;
none lists both triangles, so there is no file with a different convention (scan_convention()).

Pass rules (a failure is never dropped from the ledger):
  pass  status optimal, |obj - ref| <= 1e-6 * max(1, |ref|), row and bound violation <= 1e-6
        (violation measured by the independent recheck, relative to 1 + |bound|; absolute values also recorded)
  FAIL  anything else (other status, wrong objective, infeasible point, crash, no HiGHS reference)

Usage: python3 benchmarks/qp_maros_meszaros.py --engine ./taral --time-limit 60 [--only NAME ...]
Writes results/qp_maros_meszaros/ledger.csv (and caches downloads under out/maros_meszaros/).
"""
import argparse, csv, hashlib, json, os, subprocess, sys, time, urllib.request

import numpy as np
import highspy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ipm_check import highs_solve, rel_err  # shared helpers: HiGHS solve of the same file, |a-b|/max(1,|b|)

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR = os.path.join(ROOT, "out", "maros_meszaros")
SOURCES = os.path.join(ROOT, "results", "qp_maros_meszaros", "sources.csv")
OUT = os.path.join(ROOT, "results", "qp_maros_meszaros", "ledger.csv")
INF = 1e20
OBJ_TOL = 1e-6
VIOL_TOL = 1e-6


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def instances(max_size):
    """(name, expected sha256, url) for every .QPS file in sources.csv with rows+cols <= max_size."""
    out = []
    for r in csv.DictReader(open(SOURCES)):
        if r["file"].endswith(".QPS") and int(r["rows"]) + int(r["cols"]) <= max_size:
            out.append((r["file"][:-4], r["sha256"], r["url"]))
    return out


def fetch(name, want, url):
    path = os.path.join(DIR, name + ".mps")  # HiGHS picks the reader by extension
    os.makedirs(DIR, exist_ok=True)
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    if sha(path) != want:
        raise SystemExit("sha256 mismatch for %s: %s" % (name, sha(path)))
    return path


def scan_convention(path):
    """Files listing both (i,j) and (j,i) off-diagonals in QUADOBJ would differ in convention; returns that count."""
    sec, seen, both = None, set(), 0
    for line in open(path, errors="replace"):
        if line[0] not in " \t":
            sec = line.split()[0] if line.strip() else sec
            continue
        t = line.split()
        if sec == "QUADOBJ" and len(t) >= 3 and t[0] != t[1]:
            both += (t[1], t[0]) in seen
            seen.add((t[0], t[1]))
    return both


def highs_model(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(path)
    return h


def check_point(path, sol):
    """Independent recheck on the original model read by HiGHS: row/bound violations (absolute and relative to
    1 + |bound|) and the objective with the Hessian term. Returns None when the .sol does not cover the columns."""
    m = highs_model(path).getModel()
    lp = m.lp_
    vals = {}
    for line in open(sol):
        t = line.rsplit(None, 1)
        if len(t) == 2:
            vals[t[0].strip()] = float(t[1])
    if any(c not in vals for c in lp.col_names_):
        return None
    x = np.array([vals[c] for c in lp.col_names_])
    a = lp.a_matrix_
    start, index, value = np.array(a.start_), np.array(a.index_), np.array(a.value_)
    act = np.zeros(lp.num_row_)
    for j in range(lp.num_col_):
        s, e = start[j], start[j + 1]
        act[index[s:e]] += value[s:e] * x[j]

    def viol(v, lo, up):
        ab = rl = 0.0
        for vi, l, u in zip(v, lo, up):
            for d, b in ((l - vi, l) if l > -INF else (0, 0), (vi - u, u) if u < INF else (0, 0)):
                ab, rl = max(ab, d), max(rl, d / (1 + abs(b)))
        return ab, rl

    row_abs, row_rel = viol(act, lp.row_lower_, lp.row_upper_)
    bnd_abs, bnd_rel = viol(x, lp.col_lower_, lp.col_upper_)
    obj = lp.offset_ + float(np.dot(lp.col_cost_, x))
    hs = m.hessian_
    if hs.dim_ > 0:  # lower triangle of symmetric Q, objective 0.5 x'Qx
        hst, hi, hv = np.array(hs.start_), np.array(hs.index_), np.array(hs.value_)
        q = 0.0
        for j in range(hs.dim_):
            for p in range(hst[j], hst[j + 1]):
                q += hv[p] * x[hi[p]] * x[j] * (1.0 if hi[p] == j else 2.0)
        obj += 0.5 * q
    return {"row_abs": row_abs, "row_rel": row_rel, "bnd_abs": bnd_abs, "bnd_rel": bnd_rel, "obj": obj,
            "rows": lp.num_row_, "cols": lp.num_col_}


def highs_ref(items, limit):
    cache_path = os.path.join(DIR, "highs_ref.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    for n, path in items:
        if n not in cache:
            r = highs_solve(path, limit)
            r["highs_version"] = highspy.Highs().version()
            cache[n] = r
            json.dump(cache, open(cache_path, "w"), indent=1)
    return cache


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=os.path.join(ROOT, "out", "taral"))
    ap.add_argument("--time-limit", type=float, default=60)
    ap.add_argument("--highs-limit", type=float, default=300)
    ap.add_argument("--max-size", type=int, default=5000, help="run instances with rows + cols <= this")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    items = [(n, fetch(n, s, u)) for n, s, u in instances(a.max_size) if not a.only or n in a.only]
    ref = highs_ref(items, a.highs_limit)
    cols = ["instance", "rows", "cols", "status", "objective", "ref_objective", "ref_status", "rel_err", "primal_res",
            "dual_res", "gap", "iterations", "chk_objective", "chk_row_viol", "chk_bound_viol", "chk_row_viol_abs",
            "chk_bound_viol_abs", "wall_s", "ref_wall_s", "verdict", "reason"]
    rows = []
    for n, path in items:
        js, sol = os.path.join(DIR, n + ".json"), os.path.join(DIR, n + ".sol")
        for p in (js, sol):
            if os.path.exists(p):
                os.remove(p)
        t0 = time.time()
        try:
            subprocess.run([a.engine, path, "--time-limit", str(a.time_limit), "--json", js, "--sol", sol],
                           capture_output=True, timeout=a.time_limit * 3 + 60)
            r = json.load(open(js))
        except Exception as e:  # crash or hang is a failure, kept in the ledger
            r = {"status": "crash", "objective": None, "wall_s": time.time() - t0, "message": repr(e)}
        hr = ref[n]
        refv = hr["objective"]
        chk = check_point(path, sol) if os.path.exists(sol) else None
        if chk is None:
            lp = highs_model(path).getLp()
            chk = {"rows": lp.num_row_, "cols": lp.num_col_}
        st, reason, verdict = r["status"], "", "FAIL"
        if scan_convention(path):
            reason += "QUADOBJ lists both triangles (different convention); "
        err = rel_err(r.get("objective"), refv)
        if refv is None:
            reason += "no HiGHS reference (%s); " % hr["status"]
        if st != "optimal":
            reason += "status %s: %s; " % (st, r.get("message", ""))
        elif "obj" not in chk:
            reason += "optimal claimed but no usable .sol; "
        else:
            if abs(chk["obj"] - r["objective"]) > OBJ_TOL * max(1, abs(chk["obj"])):
                reason += "reported objective %r != recomputed %r; " % (r["objective"], chk["obj"])
            if err is not None and err > OBJ_TOL:
                reason += "objective %r vs HiGHS %r (rel err %.3g); " % (r["objective"], refv, err)
            if chk["row_rel"] > VIOL_TOL or chk["bnd_rel"] > VIOL_TOL:
                reason += "violation row %.3g bound %.3g; " % (chk["row_rel"], chk["bnd_rel"])
            if not reason and err is not None:
                verdict = "pass"
        rows.append([n, chk["rows"], chk["cols"], st, r.get("objective"), refv, hr["status"], err,
                     r.get("primal_res"), r.get("dual_res"), r.get("gap"), r.get("iterations"), chk.get("obj"),
                     chk.get("row_rel"), chk.get("bnd_rel"), chk.get("row_abs"), chk.get("bnd_abs"),
                     round(r["wall_s"], 3), round(hr["wall_s"], 3), verdict, reason.strip("; ")])
        print(" ".join(str(v) for v in rows[-1]), flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
    tally = {}
    for row in rows:
        tally.setdefault(row[-2], []).append(row[0])
    for k, v in sorted(tally.items()):
        print("%s %d: %s" % (k, len(v), " ".join(v)))
    print("ledger:", a.out)


if __name__ == "__main__":
    main()
