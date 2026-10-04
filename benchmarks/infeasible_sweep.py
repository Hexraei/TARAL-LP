#!/usr/bin/env python3
"""Sweep of `taral --explain-infeasible` over infeasible LP models, with the independent checker.

Usage: infeasible_sweep.py --engine BIN [--corpus DIR_WITH_MPS] [--generated N] [--seed S] [--time-limit SEC]
                           [--highs-iis] [--out results.json]
No corpus is stored in the repository: point --corpus at a directory of .mps files (for example the 29 Netlib
infeasible models). Per model it records, as MEASURED values: the engine's plain status, the explanation status,
rows explained, unproven rows, wall time, whether the independent checker (infeasibility_explanation_check) accepts
the explanation, the HiGHS status, and optionally the HiGHS IIS size.

Guarantees, stated per procedure (never labelled "smallest"):
  taral explanation   row-irreducible under retained column bounds, floating point, deletion order dependent.
  HiGHS IIS           whatever HiGHS' IIS routine returns for the options used here; row/column counts are
                      reference sizes only: no claim of irreducibility or minimality is made from them.
  planted (generated) the planted conflict size k is an UPPER bound on the smallest conflict; the explained set
                      can be smaller or different, and a size equal to k proves nothing about minimality.
A model for which the engine does not return `infeasible` is reported with its status and is not explained."""
import argparse, json, math, os, random, subprocess, sys, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import infeasibility_explanation_check as chk

def run(cmd, timeout):
    t = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p, time.time() - t
    except subprocess.TimeoutExpired:
        return None, time.time() - t

def status_line(p):
    if p is None: return "harness_timeout"
    lines = [l for l in p.stdout.splitlines() if l.startswith("status ")]
    return lines[-1].split()[1] if lines else "no_status_line"

IIS_STRATEGY = 6  # kIisStrategyFromLp (2) | kIisStrategyIrreducible (4), HiGHS bit flags; recorded in the output

def highs_info(path, want_iis):
    import highspy
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.readModel(chk.normalize_mps(path))
    if want_iis: h.setOptionValue("iis_strategy", IIS_STRATEGY); h.setOptionValue("time_limit", 120.0)
    h.run()
    st = h.modelStatusToString(h.getModelStatus())
    out = {"highs_status": st}
    if want_iis and st == "Infeasible":
        try:
            stt, iis = h.getIis()
            out.update(highs_iis_rows=len(iis.row_index_), highs_iis_cols=len(iis.col_index_), highs_iis_valid=bool(iis.valid_), highs_iis_call_status=str(stt), highs_iis_strategy=IIS_STRATEGY)
        except Exception as ex:
            out["highs_iis_error"] = "%s: %s" % (type(ex).__name__, ex)
    return out

def generate(rng, path):
    """Planted-conflict LP, neutral integer seeds: a feasible random system plus one row that is a linear
    combination of k-1 earlier rows with the right-hand side shifted past feasibility. Returns k."""
    n = rng.randint(8, 30); m = rng.randint(6, 20); k = rng.randint(2, min(5, m))
    x0 = [rng.randint(0, 5) for _ in range(n)]
    rows = []
    for _ in range(m):
        idx = rng.sample(range(n), rng.randint(2, min(6, n)))
        a = {j: rng.randint(-4, 4) or 1 for j in idx}
        act = sum(v * x0[j] for j, v in a.items()); rows.append(("L", a, act + rng.randint(0, 3)))
    base = rng.sample(range(m), k - 1); comb = {}
    for r in base:
        for j, v in rows[r][1].items(): comb[j] = comb.get(j, 0) + v
    comb = {j: v for j, v in comb.items() if v}
    if not comb: comb = {0: 1}; base = [0]; k = 2
    rhs_sum = sum(rows[r][2] for r in base)
    rows.append(("G", comb, rhs_sum + rng.randint(1, 3)))  # sum(base rows) <= rhs_sum but combined row >= rhs_sum + d
    L = ["NAME GEN", "ROWS", " N OBJ"] + [" %s R%d" % (t, i) for i, (t, a, b) in enumerate(rows)] + ["COLUMNS"]
    for j in range(n):
        ents = [(i, a[j]) for i, (t, a, b) in enumerate(rows) if j in a] or [(0, 0)]
        for i, v in ents: L.append(" X%d R%d %d" % (j, i, v))
    L += ["RHS"] + [" RHS R%d %d" % (i, b) for i, (t, a, b) in enumerate(rows)] + ["BOUNDS"] + [" UP BND X%d 1000" % j for j in range(n)] + ["ENDATA"]
    open(path, "w").write("\n".join(L) + "\n")
    return len(base) + 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True); ap.add_argument("--corpus"); ap.add_argument("--generated", type=int, default=0)
    ap.add_argument("--seed", type=int, default=7); ap.add_argument("--time-limit", type=float, default=60.0)
    ap.add_argument("--highs-iis", action="store_true"); ap.add_argument("--out"); ap.add_argument("--keep", help="directory to keep per-model engine stdout and explanation JSON")
    a = ap.parse_args()
    tmp = a.keep or tempfile.mkdtemp(); os.makedirs(tmp, exist_ok=True); items = []
    if a.corpus:
        for f in sorted(os.listdir(a.corpus)):
            if f.endswith(".mps"): items.append((f[:-4], os.path.join(a.corpus, f), None))
    rng = random.Random(a.seed)
    for g in range(a.generated):
        p = os.path.join(tmp, "gen%d.mps" % g); items.append(("gen%d" % g, p, generate(rng, p)))
    results = []
    for name, path, planted in items:
        r = {"model": name, "planted_k_upper_bound": planted}
        p, w = run([a.engine, path, "--time-limit", str(a.time_limit)], a.time_limit + 30)
        r["engine_status"] = status_line(p); r["engine_wall_s"] = round(w, 3)
        if a.keep and p is not None: open(os.path.join(tmp, name + ".plain.stdout.txt"), "w").write(p.stdout)
        try: r.update(highs_info(path, a.highs_iis))
        except Exception as ex: r["highs_error"] = "%s: %s" % (type(ex).__name__, ex)
        if r["engine_status"] == "infeasible":
            ej = os.path.join(tmp, name + ".explain.json")
            p, w = run([a.engine, path, "--explain-infeasible", ej, "--time-limit", str(a.time_limit)], a.time_limit + 30)
            if a.keep and p is not None: open(os.path.join(tmp, name + ".explain.stdout.txt"), "w").write(p.stdout)
            r["explain_exit"] = None if p is None else p.returncode; r["explain_wall_s"] = round(w, 3)
            if os.path.exists(ej):
                e = json.load(open(ej)); r["explain_status"] = e.get("status"); r["rows"] = len(e.get("rows") or []); up = e.get("unproven_rows"); r["unproven_rows"] = len(up) if isinstance(up, list) else up
                bc = e.get("bound_columns"); r["bound_columns"] = len(bc) if isinstance(bc, list) else bc; r["relaxation_status"] = (e.get("relaxation") or {}).get("status")
                try:
                    st, k, errs = chk.check(path, ej)
                    relax = [x for x in errs if x.startswith(("relaxation", "objective"))]; rowerr = [x for x in errs if x not in relax]
                    r["checker"] = "PASS" if not errs else "FAIL"; r["checker_rows"] = "PASS" if not rowerr else "FAIL"; r["checker_relaxation"] = "PASS" if not relax else "FAIL"; r["checker_errors"] = errs[:3]
                except Exception as ex:
                    r["checker"] = "ERROR"; r["checker_errors"] = ["%s: %s" % (type(ex).__name__, ex)]
            else: r["explain_status"] = "no_output"
        else: r["explain_status"] = "not_run_engine_not_infeasible"
        results.append(r)
        print("%-14s engine=%-18s highs=%-12s explain=%-22s rows=%s checker=%s" % (name, r["engine_status"], r.get("highs_status"), r.get("explain_status"), r.get("rows"), r.get("checker")), flush=True)
    summ = {"models": len(results), "engine_infeasible": sum(r["engine_status"] == "infeasible" for r in results),
            "engine_other": sorted({r["model"]: r["engine_status"] for r in results if r["engine_status"] != "infeasible"}.items()),
            "explained_irreducible": sum(r.get("explain_status") == "irreducible" for r in results),
            "checker_pass": sum(r.get("checker") == "PASS" for r in results), "checker_not_pass": [r["model"] for r in results if r.get("checker") in ("FAIL", "ERROR")],
            "row_proof_not_pass": [r["model"] for r in results if r.get("checker_rows") == "FAIL" or r.get("checker") == "ERROR"],
            "relaxation_not_pass": [r["model"] for r in results if r.get("checker_relaxation") == "FAIL"]}
    print(json.dumps(summ))
    if a.out: json.dump({"summary": summ, "results": results, "engine": a.engine, "seed": a.seed, "time_limit": a.time_limit}, open(a.out, "w"), indent=1)
    return 1 if summ["row_proof_not_pass"] else 0  # relaxation-objective disagreements are reported, not hidden, but do not hide a verified row proof
sys.exit(main())
