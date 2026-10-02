"""Shared tooling for the IPM tests: run the IPM driver, re-check its solution independently, solve with HiGHS.

recheck() reads the ORIGINAL MPS file with HiGHS's reader (not TARAL's parser), maps the IPM's .sol
vector by column name and recomputes row activities, bound/row violations and the objective
(c'x + 0.5 x'Qx + offset) from scratch. Build the driver first:
  g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o out/ipm_run out/ipm_main.cpp \
      src/ipm.cpp src/mps.cpp
"""
import json, subprocess, time
from pathlib import Path
import numpy as np
import highspy

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "out" / "ipm_run"
INF = 1e20


def run_ipm(mps, tag_dir, time_limit=60.0, extra=()):
    tag_dir = Path(tag_dir)
    tag_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(mps).name.replace(".mps", "").replace(".qps", "").replace(".SIF", "")
    js, sol = tag_dir / f"{stem}.json", tag_dir / f"{stem}.sol"
    for f in (js, sol):
        f.unlink(missing_ok=True)
    t = time.time()
    try:
        subprocess.run([str(BIN), str(mps), "--time-limit", str(time_limit), "--json", str(js), "--sol", str(sol), *extra],
                       capture_output=True, timeout=time_limit + 60)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "wall_s": time.time() - t}, None
    if not js.exists():
        return {"status": "crash", "wall_s": time.time() - t}, None
    return json.loads(js.read_text()), (sol if sol.exists() else None)


def _model(mps):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(str(mps))
    return h


def recheck(mps, sol):
    """Independent primal recheck of the .sol x on the original model. Returns dict or None."""
    h = _model(mps)
    m = h.getModel()
    lp = m.lp_
    n, mrows = lp.num_col_, lp.num_row_
    names = list(lp.col_names_)
    # Columns are matched by position (both readers keep first-appearance order); names are
    # compared too, ignoring whitespace, since fixed-column names may contain spaces.
    xs, sn = [], []
    for line in Path(sol).read_text().splitlines():
        t = line.split()
        if t and t[0] == "C":
            xs.append(float(t[-2])); sn.append("".join(t[1:-2]))
    if len(xs) != n or any(a != "".join(b.split()) for a, b in zip(sn, names)):
        raise ValueError(f"column mismatch between readers in {mps}")
    x = np.array(xs)
    a = lp.a_matrix_
    start, idx, val = np.array(a.start_), np.array(a.index_), np.array(a.value_)
    act = np.zeros(mrows)
    for j in range(n):
        s, e = start[j], start[j + 1]
        act[idx[s:e]] += val[s:e] * x[j]
    rl, ru = np.array(lp.row_lower_), np.array(lp.row_upper_)
    cl, cu = np.array(lp.col_lower_), np.array(lp.col_upper_)
    def viol(v, lo, up):
        lo = np.where(lo <= -INF, -np.inf, lo); up = np.where(up >= INF, np.inf, up)
        a = np.maximum(np.maximum(lo - v, v - up), 0.0)
        bnd = np.where(v < lo, lo, up)
        rel = np.where(a > 0, a / (1 + np.abs(np.where(np.isfinite(bnd), bnd, 0))), 0)
        return (float(a.max()) if a.size else 0.0), (float(rel.max()) if rel.size else 0.0)
    rv, rrel = viol(act, rl, ru)
    bv, brel = viol(x, cl, cu)
    obj = float(np.dot(lp.col_cost_, x)) + lp.offset_
    hs = m.hessian_
    if hs.dim_ > 0:
        hst, hi, hv = np.array(hs.start_), np.array(hs.index_), np.array(hs.value_)
        q = 0.0
        for j in range(hs.dim_):
            for p in range(hst[j], hst[j + 1]):
                i = hi[p]
                q += hv[p] * x[i] * x[j] * (1.0 if i == j else 2.0)  # lower triangle of symmetric Q
        obj += 0.5 * q
    return {"chk_obj": obj, "chk_row_viol": rv, "chk_row_rel": rrel, "chk_bound_viol": bv, "chk_bound_rel": brel}


def highs_solve(mps, time_limit=60.0):
    h = _model(mps)
    h.setOptionValue("time_limit", float(time_limit))
    t = time.time()
    try:
        h.run()
    except Exception as e:  # HiGHS raises on e.g. nonconvex QP in some versions
        return {"status": "error:" + str(e)[:60], "objective": None, "wall_s": time.time() - t}
    st = h.modelStatusToString(h.getModelStatus())
    obj = h.getInfo().objective_function_value if h.getModelStatus() == highspy.HighsModelStatus.kOptimal else None
    return {"status": st, "objective": obj, "wall_s": time.time() - t}


def rel_err(a, b):
    if a is None or b is None:
        return None
    return abs(a - b) / max(1.0, abs(b))
