#!/usr/bin/env python3
"""Differential harness: engine binary vs HiGHS on every generated case, with an independent exact checker.

Per case: write the MPS, run the engine (time cap per case), parse the ORIGINAL file with the oracle parser,
solve the same model with HiGHS (presolve on; on any disagreement re-solve with presolve off before blaming
the engine), and check the engine's returned point against the original model in exact arithmetic.

Gates (all must hold for PASS):
  status equal (optimal / infeasible / unbounded)
  |obj_ours - obj_ref| <= 1e-6 * max(1, |obj_ref|)   (obj_ours recomputed exactly from the returned point)
  row and bound violation of the returned point <= 1e-6 * (1 + |violated bound|)  (absolute value also recorded)
  MILP: integrality violation <= 1e-6
Usage: harness.py --engine ./taral --out DIR [--cats a,b] [--limit N] [--jobs J] [--time-limit 60]
                  [--kind lp|milp|all|big|biglp|lp+big|lp+biglp|all+big] [--method simplex|dual|ipm]
--method is passed to the engine as `--method M` (simplex = no flag, the default path). The engine only honours it for
pure LPs (MILPs are routed to branch and bound whatever the flag says), so dual/ipm runs normally use --kind lp.
"""
import argparse
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from mpsio import parse_mps, exact_check, read_sol, OracleError, INF
from ref import solve_ref

TOL_OBJ = 1e-6
TOL_VIOL = 1e-6
DETERMINATE = ("optimal", "infeasible", "unbounded")


def run_engine(engine, path, sol, js, tl, method="simplex"):
    t0 = time.time()
    cmd = [engine, path, "--time-limit", str(tl), "--sol", sol, "--json", js]
    if method != "simplex":
        cmd += ["--method", method]
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=tl + 30)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        return dict(status="harness_timeout", objective=None, wall=time.time() - t0, rc=None, message="killed")
    wall = time.time() - t0
    try:
        j = json.load(open(js))
    except Exception:
        return dict(status="no_json_rc%s" % rc, objective=None, wall=wall, rc=rc, message=p.stderr.decode()[:200])
    j["wall"] = wall
    j["rc"] = rc
    return j


def norm_status(s):
    if s in ("unbounded_relaxation",):
        return "unbounded_relaxation"
    return s


def rel_ok(a, b):
    return abs(a - b) <= TOL_OBJ * max(1.0, abs(b))


def reference(m, tl, expect):
    """Return (status, obj, x, notes, ref1). Presolve on first; ambiguous/odd statuses are re-solved presolve off."""
    notes = []
    r1 = solve_ref(m, presolve=True, time_limit=tl)
    st = r1["status"]
    if st in DETERMINATE or st in ("ref_hang", "ref_crash"):
        return r1, notes, r1
    r2 = solve_ref(m, presolve=False, time_limit=tl)
    notes.append("ref presolve-on status %s -> presolve-off %s" % (st, r2["status"]))
    return r2, notes, r1


def judge(case, m, ours, sol_path, tl):
    rec = dict(ours_status=ours.get("status"), ours_obj=ours.get("objective"), ours_wall=round(ours.get("wall", 0), 4))
    notes = []
    verdict = "PASS"
    ostat = norm_status(ours.get("status"))
    if ostat == "parse_error":
        rec["verdict"], rec["detail"] = "PARSE_ERR", ours.get("message", "")
        return rec
    t0 = time.time()
    ref, rnotes, ref1 = reference(m, tl, case.expect)
    notes += rnotes
    rec["ref_status"], rec["ref_obj"] = ref["status"], ref["objective"]
    rec["ref_wall"] = round(time.time() - t0, 4)
    rstat = ref["status"]
    if case.expect and rstat in DETERMINATE and rstat != case.expect:
        # The construction proves the status; the reference disagrees with it. Score the engine against the construction.
        r_off = solve_ref(m, presolve=False, time_limit=tl) if ref is ref1 else ref
        if r_off["status"] == case.expect:
            ref, rstat = r_off, r_off["status"]
            rec["ref_status"], rec["ref_obj"] = ref["status"], ref["objective"]
            notes.append("reference presolve gave %s; presolve-off agrees with the construction" % rec["ref_status"])
        else:
            rec["construction_vs_ref"] = "constructed %s, reference %s (presolve-off %s)" % (case.expect, rstat, r_off["status"])
            if ostat == case.expect:
                if ostat == "optimal":
                    try:
                        chk0 = exact_check(m, read_sol(sol_path, m))
                    except Exception:
                        chk0 = None
                    if chk0 is None or chk0["viol_row_rel"] > TOL_VIOL or chk0["viol_bnd_rel"] > TOL_VIOL or chk0["int_viol"] > TOL_VIOL:
                        rec["verdict"], rec["detail"] = "FAIL_FEAS", "engine optimal but point infeasible; " + rec["construction_vs_ref"]
                        return rec
                    rec["ours_obj_recomputed"] = chk0["objective"]
                rec["verdict"], rec["detail"] = "PASS_REF", "engine matches the construction; " + rec["construction_vs_ref"]
                return rec
            rec["verdict"] = "FAIL_STATUS" if ostat in DETERMINATE else "NOANS"
            rec["detail"] = "engine %s (%s); %s" % (ostat, ours.get("message", "")[:80], rec["construction_vs_ref"])
            return rec
    chk = None
    if ostat == "optimal":
        try:
            x = read_sol(sol_path, m)
            chk = exact_check(m, x)
        except Exception as e:  # missing/garbled .sol
            rec["verdict"], rec["detail"] = "FAIL_NOSOL", repr(e)
            return rec
        rec["viol_row_abs"], rec["viol_row_rel"] = chk["viol_row_abs"], chk["viol_row_rel"]
        rec["viol_bnd_abs"], rec["viol_bnd_rel"] = chk["viol_bnd_abs"], chk["viol_bnd_rel"]
        rec["int_viol"] = chk["int_viol"]
        rec["ours_obj_recomputed"] = chk["objective"]
    # ---- status comparison
    if rstat in ("ref_hang", "ref_crash"):
        rec["verdict"], rec["detail"] = "REF_UNRESOLVED", "reference solver %s (HiGHS 1.15.1); engine said %s%s" % (
            rstat, ostat, " obj %r, exact checks %s" % (chk["objective"], "pass" if chk["viol_row_rel"] <= TOL_VIOL and chk["viol_bnd_rel"] <= TOL_VIOL and chk["int_viol"] <= TOL_VIOL else "FAIL") if chk else "")
        return rec
    if rstat not in DETERMINATE:
        if ostat in DETERMINATE and case.expect == ostat:
            rec["verdict"], rec["detail"] = "PASS_REF", "reference %s; engine matches construction" % rstat
            return rec
        rec["verdict"], rec["detail"] = "REF_UNRESOLVED", "reference status %s, engine %s" % (rstat, ostat)
        return rec
    if ostat != rstat:
        # arbitration
        r_off = solve_ref(m, presolve=False, time_limit=tl) if ref1["status"] == rstat and ref is ref1 else ref
        if r_off["status"] == ostat and ostat in DETERMINATE:
            rec["verdict"], rec["detail"] = "PASS_REF", "reference presolve gave %s; presolve-off agrees with engine (%s)" % (rstat, ostat)
            if ostat == "optimal":
                if not (rel_ok(chk["objective"], r_off["objective"])):
                    rec["verdict"], rec["detail"] = "FAIL_OBJ", "status ok vs presolve-off but objective %r vs %r" % (chk["objective"], r_off["objective"])
            return rec
        if ostat == "optimal" and chk and chk["viol_row_rel"] <= TOL_VIOL and chk["viol_bnd_rel"] <= TOL_VIOL \
                and chk["int_viol"] <= TOL_VIOL and rstat == "infeasible":
            rec["verdict"], rec["detail"] = "PASS_REF", "engine point is feasible on the original model (exact check); reference said infeasible"
            return rec
        if rstat == "optimal" and ostat == "infeasible":
            rx = ref["x"]
            rchk = exact_check(m, rx) if rx is not None else None
            if rchk and rchk["viol_row_rel"] <= TOL_VIOL and rchk["viol_bnd_rel"] <= TOL_VIOL and rchk["int_viol"] <= TOL_VIOL:
                rec["verdict"] = "FAIL_STATUS"
                rec["detail"] = "engine infeasible; reference optimal %r with a point that is feasible on the original model (exact check)" % ref["objective"]
                return rec
        kind = "NOANS" if ostat not in DETERMINATE else "FAIL_STATUS"
        rec["verdict"], rec["detail"] = kind, "engine %s (%s) vs reference %s (presolve-off %s)" % (ostat, ours.get("message", "")[:80], rstat, r_off["status"])
        return rec
    # ---- statuses agree
    if ostat != "optimal":
        rec["verdict"], rec["detail"] = "PASS", ostat
        return rec
    fails = []
    if chk["viol_row_rel"] > TOL_VIOL:
        fails.append("row violation %.3g (abs %.3g)" % (chk["viol_row_rel"], chk["viol_row_abs"]))
    if chk["viol_bnd_rel"] > TOL_VIOL:
        fails.append("bound violation %.3g (abs %.3g)" % (chk["viol_bnd_rel"], chk["viol_bnd_abs"]))
    if case.mip and chk["int_viol"] > TOL_VIOL:
        fails.append("integrality violation %.3g" % chk["int_viol"])
    if fails:
        rec["verdict"], rec["detail"] = "FAIL_FEAS", "; ".join(fails)
        return rec
    if ours.get("objective") is not None and not rel_ok(ours["objective"], chk["objective"]):
        rec["verdict"], rec["detail"] = "FAIL_OBJ", "reported objective %r differs from recomputed %r" % (ours["objective"], chk["objective"])
        return rec
    if not rel_ok(chk["objective"], ref["objective"]):
        better = (chk["objective"] < ref["objective"]) != m.maximize
        if case.exact_obj is not None and rel_ok(chk["objective"], case.exact_obj):
            rec["verdict"], rec["detail"] = "PASS_REF", "engine matches the exact optimum %r; reference %r is inexact" % (case.exact_obj, ref["objective"])
            return rec
        r_off = solve_ref(m, presolve=False, time_limit=tl)
        if r_off["status"] == "optimal" and rel_ok(chk["objective"], r_off["objective"]):
            rec["verdict"], rec["detail"] = "PASS_REF", "reference presolve objective %r; presolve-off agrees with engine" % ref["objective"]
            return rec
        if better:
            rec["verdict"] = "PASS_REF"
            rec["detail"] = "engine objective %r is better than reference %r and its point is feasible (exact check)" % (chk["objective"], ref["objective"])
            return rec
        rec["verdict"], rec["detail"] = "FAIL_OBJ", "engine %r vs reference %r (rel %.3g)" % (chk["objective"], ref["objective"], abs(chk["objective"] - ref["objective"]) / max(1, abs(ref["objective"])))
        if case.exact_obj is not None:
            rec["detail"] += "; exact optimum %r (engine off by %.3g, reference off by %.3g)" % (
                case.exact_obj, abs(chk["objective"] - case.exact_obj) / max(1, abs(case.exact_obj)),
                abs(ref["objective"] - case.exact_obj) / max(1, abs(case.exact_obj)))
        return rec
    rec["verdict"] = "PASS"
    rec["detail"] = "; ".join(notes)
    return rec


def work(a):
    engine, outdir, cat, i, tl, method = a
    case = gen.gen(cat, i)
    d = os.path.join(outdir, "cases", cat)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, case.id + ".mps")
    with open(path, "w") as f:
        f.write(case.text)
    sol, js = path[:-4] + ".sol", path[:-4] + ".json"
    ours = run_engine(engine, path, sol, js, tl, method)
    rec = dict(id=case.id, cat=cat, seed=case.seed, mip=case.mip, expect=case.expect, method=method)
    try:
        m = parse_mps(path)
    except OracleError as e:
        rec.update(verdict="ORACLE_ERR", detail=str(e), ours_status=ours.get("status"))
        return rec
    try:
        rec.update(judge(case, m, ours, sol, tl))
    except Exception as e:  # harness bug must never look like a pass
        import traceback
        rec.update(verdict="HARNESS_ERR", detail=traceback.format_exc()[-400:])
    if rec["verdict"] in ("PASS", "PASS_REF"):  # keep disk small: only failures keep their files
        for p in (path, sol, js):
            if os.path.exists(p):
                os.remove(p)
    return rec


def worker_main():
    """Long-lived worker: reads 'cat index' lines on stdin, writes one JSON result line per case on stdout."""
    engine, outdir, tl = sys.argv[2], sys.argv[3], float(sys.argv[4])
    method = sys.argv[5] if len(sys.argv) > 5 else "simplex"
    for line in sys.stdin:
        cat, i = line.split()
        rec = work((engine, outdir, cat, int(i), tl, method))
        sys.stdout.write(json.dumps(rec, sort_keys=True) + "\n")
        sys.stdout.flush()


def supervise(jobs, nworkers, engine, outdir, tl, results, method="simplex"):
    """Own worker pool: a worker that dies or hangs (e.g. a crash inside the reference solver) is replaced and the
    case it held is recorded as HARNESS_ERR instead of hanging the whole run."""
    import queue
    import select
    import threading
    q = queue.Queue()
    for j in jobs:
        q.put(j)
    lock = threading.Lock()
    done = [0]
    t0 = time.time()
    cap = 2 * tl + 120

    def spawn():
        return subprocess.Popen([sys.executable, os.path.abspath(__file__), "--worker", engine, outdir, str(tl), method],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)

    def loop():
        w = spawn()
        while True:
            try:
                _, _, cat, i, _, _ = q.get_nowait()
            except queue.Empty:
                break
            rec = None
            try:
                w.stdin.write("%s %d\n" % (cat, i))
                w.stdin.flush()
                r, _, _ = select.select([w.stdout], [], [], cap)
                line = w.stdout.readline() if r else ""
                if line:
                    rec = json.loads(line)
            except Exception:
                rec = None
            if rec is None:
                seed = gen.BASE[cat] + i
                rec = dict(id="%s-%d" % (cat, seed), cat=cat, seed=seed, mip=cat in gen.IS_MIP, verdict="HARNESS_ERR", method=method,
                           detail="worker died or hung on this case (crash/hang in engine wrapper or reference solver)")
                try:
                    w.kill()
                except Exception:
                    pass
                w = spawn()
            with lock:
                results.append(rec)
                done[0] += 1
                if done[0] % 250 == 0:
                    print("  ... %d/%d cases, %.0fs" % (done[0], len(jobs), time.time() - t0), file=sys.stderr, flush=True)
        try:
            w.stdin.close()
            w.wait(timeout=10)
        except Exception:
            w.kill()

    ts = [threading.Thread(target=loop) for _ in range(nworkers)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cats", default="")
    ap.add_argument("--kind", default="all", choices=["all", "lp", "milp", "big", "lp+big", "biglp", "lp+biglp", "all+big"])
    ap.add_argument("--only", default="", help="comma-separated case ids (cat-seed) to run instead of whole categories")
    ap.add_argument("--method", default="simplex", choices=["simplex", "dual", "ipm"])
    ap.add_argument("--limit", type=int, default=0, help="cases per category (0 = full count)")
    ap.add_argument("--jobs", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    ap.add_argument("--time-limit", type=float, default=60.0)
    a = ap.parse_args()
    cats = [c for c, _, _ in gen.ALL]
    if a.kind == "lp":
        cats = [c for c, _, _ in gen.LP_CATS]
    elif a.kind == "milp":
        cats = [c for c, _, _ in gen.MILP_CATS]
    elif a.kind == "big":
        cats = [c for c, _, _ in gen.BIG_CATS]
    elif a.kind == "lp+big":
        cats = [c for c, _, _ in gen.LP_CATS + gen.BIG_CATS]
    elif a.kind == "biglp":
        cats = [c for c, _, _ in gen.BIG_CATS if not gen.IS_MIP.get(c)]
    elif a.kind == "lp+biglp":
        cats = [c for c, _, _ in gen.LP_CATS + gen.BIG_CATS if not gen.IS_MIP.get(c)]
    elif a.kind == "all+big":
        cats = [c for c, _, _ in gen.ALL + gen.BIG_CATS]
    if a.cats:
        cats = [c for c in a.cats.split(",") if c]
    jobs = []
    for c in cats:
        cnt = gen.COUNT[c] if not a.limit else min(a.limit, gen.COUNT[c])
        jobs += [(os.path.abspath(a.engine), os.path.abspath(a.out), c, i, a.time_limit, a.method) for i in range(cnt)]
    if a.only:
        jobs = []
        for cid in a.only.split(","):
            cat, seed = cid.rsplit("-", 1)
            jobs.append((os.path.abspath(a.engine), os.path.abspath(a.out), cat, int(seed) - gen.BASE[cat], a.time_limit, a.method))
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    recs = []
    supervise(jobs, a.jobs, os.path.abspath(a.engine), os.path.abspath(a.out), a.time_limit, recs, a.method)
    recs.sort(key=lambda r: (r["cat"], r["seed"]))
    with open(os.path.join(a.out, "results.jsonl"), "w") as f:
        for r in recs:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    print("harness wall %.1fs for %d cases" % (time.time() - t0, len(recs)), file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker_main()
    else:
        main()
