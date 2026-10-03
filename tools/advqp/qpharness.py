#!/usr/bin/env python3
"""Adversarial QP harness: engine (QUADOBJ -> interior point) vs HiGHS QP and vs independent KKT certificates.

Per case: write the MPS, run the engine, then judge it against (1) the construction when the generator knows the
answer (a KKT-built optimum, a verified recession ray, an arithmetic infeasibility), (2) HiGHS built from the same
data with passModel, and (3) checks computed only from the returned point: exact rational row/bound violations
(absolute and relative) and the best KKT stationarity residual (qpcheck.kkt_certificate).

Gates for PASS on an optimal answer (all must hold):
  exact row/bound violation <= 1e-6 * (1 + |violated bound|)           (absolute value always recorded)
  best KKT stationarity residual <= 1e-6   (per-column relative to 1+|g_j|; absolute recorded too) and
  certified duality gap <= 1e-6 * (1 + |objective|)   (qpcheck.dual_certificate: LP over ALL multipliers of the point)
  [reported, not gated: act_rel, the stationarity residual when only constraints within 1e-6 of a bound carry
   multipliers, i.e. how far an interior-point answer sits from the optimal face]
  |objective(x) - target| <= 1e-6 * max(1, |target|), objective recomputed exactly from x; target is the constructed
  optimum when known, else HiGHS. Engine-reported objective must match the recomputed one within the same tolerance.
Outcomes: pass, pass* (engine certified, reference or construction shown inexact), partial (correct but weaker verdict:
dual_infeasible), border (Q within numerical noise of semidefinite), wrong, noans, other.

Usage: qpharness.py --engine ./taral --out DIR [--cats a,b] [--limit N] [--jobs J] [--time-limit 30] [--selftest]
"""
import argparse
import json
import math
import os
import select
import subprocess
import sys
import threading
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qpcase
import qpcheck
import qpequality

TOL_VIOL = 1e-6
TOL_OBJ = 1e-6
TOL_KKT = 1e-6
INF = math.inf
WRONG = ("FAIL_STATUS", "FAIL_OBJ", "FAIL_FEAS", "FAIL_KKT", "FAIL_NOSOL", "PARSE_ERR")
OTHER = ("REF_UNRESOLVED", "GEN_MISMATCH", "HARNESS_ERR")


def run_engine(engine, path, sol, js, tl):
    t0 = time.time()
    for p in (sol, js):
        if os.path.exists(p):
            os.remove(p)
    try:
        pr = subprocess.run([engine, path, "--time-limit", str(tl), "--sol", sol, "--json", js],
                            capture_output=True, timeout=tl + 30)
        rc = pr.returncode
    except subprocess.TimeoutExpired:
        return dict(status="harness_timeout", objective=None, wall=time.time() - t0, message="killed")
    wall = time.time() - t0
    try:
        j = json.load(open(js))
    except Exception:
        return dict(status="no_json_rc%s" % rc, objective=None, wall=wall, message=pr.stderr.decode(errors="replace")[:200])
    j["wall"] = wall
    return j


def read_sol(path, n):
    vals = {}
    if os.path.exists(path):
        for ln in open(path):
            p = ln.split()
            if len(p) == 2:
                vals[p[0]] = float(p[1])
    if not vals:
        return None
    return np.array([vals.get("C%d" % j, 0.0) for j in range(n)])


def convexity(case):
    lam, scale = case.eig()
    s = max(1.0, scale)
    if lam >= -1e-13 * s:
        return "convex", lam
    if lam <= -1e-6 * s:
        return "nonconvex", lam
    return "border", lam


def certify_point(case, x):
    pc = qpcheck.primal_check(case, x)
    kc = qpcheck.kkt_certificate(case, x)
    return pc, kc


def judge(case, ours, x, tl):
    rec = dict(ours_status=ours.get("status"), ours_obj=ours.get("objective"), ours_wall=round(ours.get("wall", 0), 4),
               self_row=ours.get("max_row_viol"), self_bnd=ours.get("max_bound_viol"),
               self_dual=ours.get("dual_res"), self_gap=ours.get("gap"), msg=(ours.get("message") or "")[:160])
    ost = ours.get("status")
    conv, lam = convexity(case)
    rec["convexity"], rec["lam_min"] = conv, lam
    equality_model = qpequality.reduced_case(case)
    feas_convex = case.cat == "indef_convex_feasible" or equality_model is not None
    if equality_model is not None:
        rec.update(equality_convex=True, equality_lam_min=equality_model[1], equality_dimension=equality_model[2])
    if ost == "parse_error":
        rec.update(verdict="PARSE_ERR", detail=rec["msg"])
        return rec
    nonconv_model = conv == "nonconvex" and not feas_convex

    # ---- construction sanity: the constructed optimum must itself certify (generator + certificate self-check)
    if case.x0 is not None and case.expect_status == "optimal" and case.cat != "psd_free_null" and case.cat != "ill_cond" \
            and not nonconv_model and conv != "border":
        pc0, kc0 = certify_point(case, case.x0)
        rec["x0_row_abs"], rec["x0_stat_rel"] = max(pc0["row_abs"], pc0["bnd_abs"]), kc0["stat_rel"]
        if max(pc0["row_rel"], pc0["bnd_rel"]) > 1e-9 or not (kc0["stat_rel"] <= 1e-9):
            rec.update(verdict="GEN_MISMATCH", detail="constructed optimum fails its own certificate: viol %.2e/%.2e stat %.2e" % (
                pc0["row_rel"], pc0["bnd_rel"], kc0["stat_rel"]))
            return rec

    # ---- nonconvex model: only an honest 'nonconvex' (or a refusal) is acceptable
    if nonconv_model:
        if ost == "nonconvex":
            rec.update(verdict="PASS", detail="")
        elif ost == "optimal" and x is not None:
            pc, kc = certify_point(case, x)
            rec.update(verdict="FAIL_STATUS", detail="optimal claimed on indefinite Q (lam_min %.3g); point is %s, stationarity %.2e" % (
                lam, "feasible" if max(pc["row_rel"], pc["bnd_rel"]) <= TOL_VIOL else "infeasible by %.2e" % max(pc["row_rel"], pc["bnd_rel"]), kc["stat_rel"]))
            rec.update(row_abs=pc["row_abs"], bnd_abs=pc["bnd_abs"], stat_rel=kc["stat_rel"], stat_abs=kc["stat_abs"])
        elif ost in ("infeasible", "unbounded", "dual_infeasible"):
            rec.update(verdict="FAIL_STATUS", detail="%s on a feasible boxed model with indefinite Q (feasible point constructed)" % ost)
        else:
            rec.update(verdict="NOANS", detail=ost)
        return rec

    # ---- reference
    reference_case = equality_model[0] if equality_model is not None else case
    ref = qpcheck.reference(reference_case, presolve=True, tl=tl)
    if ref["status"] not in ("optimal", "infeasible", "unbounded"):
        ref2 = qpcheck.reference(reference_case, presolve=False, tl=tl)
        rec["ref_note"] = "presolve-on %s -> presolve-off %s" % (ref["status"], ref2["status"])
        if ref2["status"] in ("optimal", "infeasible", "unbounded"):
            ref = ref2
    rec["ref_status"], rec["ref_obj"] = ref["status"], ref["objective"]
    truth = case.expect_status
    target = case.expect_obj if truth == "optimal" and case.expect_obj is not None else None
    if truth is None:
        truth = ref["status"] if ref["status"] in ("optimal", "infeasible", "unbounded") else None
    if truth == "optimal" and target is None and ref["status"] == "optimal":
        target = ref["objective"]
    if truth == "optimal" and ref["status"] == "optimal" and case.expect_obj is not None and not feas_convex:
        d = abs(ref["objective"] - case.expect_obj) / max(1.0, abs(case.expect_obj))
        if d > 1e-6:
            rec["ref_off_construction"] = d  # HiGHS disagrees with the constructed optimum; the engine is judged vs construction
    if case.expect_status == "optimal" and target is None and ref["status"] == "optimal":
        target = ref["objective"]
    rec["truth"], rec["target"] = truth, target
    sgn = -1.0 if case.maximize else 1.0

    if ost == "nonconvex":
        if conv == "border" and equality_model is None:
            rec.update(verdict="BORDER", detail="refused as nonconvex, lam_min %.3g" % lam)
        elif feas_convex:
            rec.update(verdict="NOANS", detail="refused: Q indefinite although convex on {Ax=b} (engine requires PSD Q); lam_min %.3g" % lam)
        else:
            rec.update(verdict="FAIL_STATUS", detail="nonconvex reported for a PSD Q (lam_min %.3g)" % lam)
        return rec

    if ost == "optimal":
        if x is None:
            rec.update(verdict="FAIL_NOSOL", detail="optimal but no solution file")
            return rec
        pc, kc = certify_point(case, x)
        vr = max(pc["row_rel"], pc["bnd_rel"])
        rec.update(row_abs=pc["row_abs"], bnd_abs=pc["bnd_abs"], row_rel=pc["row_rel"], bnd_rel=pc["bnd_rel"],
                   obj_exact=pc["objective"], stat_rel=kc["stat_rel"], stat_abs=kc["stat_abs"], gap=kc["gap"],
                   gap_rel=kc["gap"] / (1.0 + abs(pc["objective"])), act_rel=kc["act_rel"], act_abs=kc["act_abs"],
                   n_act=kc["n_act"], g_inf=kc["g_inf"])
        if case.x0 is not None and conv == "convex" and lam > 1e-9 * max(1.0, case.eig()[1]) and case.expect_obj is not None and case.x0 is not None:
            rec["x_err"] = float(np.max(np.abs(np.asarray(x) - case.x0)))  # unique minimiser (PD Q): distance to it
        rep = ours.get("objective")
        probs = []
        if rep is None or abs(rep - pc["objective"]) > TOL_OBJ * max(1.0, abs(pc["objective"])):
            probs.append(("FAIL_OBJ", "reported objective %r vs exact %r" % (rep, pc["objective"])))
        if vr > TOL_VIOL:
            probs.append(("FAIL_FEAS", "violation rel %.3g (abs row %.3g bnd %.3g)" % (vr, pc["row_abs"], pc["bnd_abs"])))
        if not (kc["stat_rel"] <= TOL_KKT):
            probs.append(("FAIL_KKT", "stationarity residual %.3g rel / %.3g abs" % (kc["stat_rel"], kc["stat_abs"])))
        elif not (rec["gap_rel"] <= TOL_KKT):
            probs.append(("FAIL_KKT", "certified duality gap %.3g (rel %.3g)" % (kc["gap"], rec["gap_rel"])))
        better = False
        if truth in ("infeasible", "unbounded"):
            if case.expect_status in ("infeasible", "unbounded"):
                probs.append(("FAIL_STATUS", "optimal on a model that is %s by construction" % truth))
            elif vr <= TOL_VIOL and kc["stat_rel"] <= TOL_KKT and rec["gap_rel"] <= TOL_KKT:
                better = True  # exactly feasible + KKT-certified point: the reference verdict is wrong
                rec["pass_ref_why"] = "reference says %s; engine point is feasible and KKT-certified (convex: global optimum)" % truth
            else:
                probs.append(("FAIL_STATUS", "optimal but reference says %s" % truth))
        elif target is not None:
            gap = sgn * (pc["objective"] - target)  # >0: worse than target
            rec["obj_gap"] = gap
            if gap > TOL_OBJ * max(1.0, abs(target)):
                probs.append(("FAIL_OBJ", "objective %.12g worse than target %.12g (rel %.2e)" % (pc["objective"], target, gap / max(1.0, abs(target)))))
            elif gap < -TOL_OBJ * max(1.0, abs(target)):
                better = True
                rec["pass_ref_why"] = "engine objective %.12g beats %s %.12g" % (pc["objective"], "construction" if case.expect_obj is not None and not feas_convex else "reference", target)
        elif truth is None:
            rec.update(verdict="REF_UNRESOLVED", detail="reference %s; engine point feasible=%s stationarity %.2e" % (ref["status"], vr <= TOL_VIOL, kc["stat_rel"]))
            if not probs:
                return rec
        if probs:
            order = [p[0] for p in probs]
            first = next(p for p in probs if p[0] in ("FAIL_STATUS", "FAIL_FEAS", "FAIL_OBJ", "FAIL_KKT"))
            rec.update(verdict=first[0], detail="; ".join(p[1] for p in probs), all_fail=order)
            return rec
        if conv == "border" and equality_model is None:
            rec.update(verdict="BORDER", detail="accepted a Q within noise of semidefinite (lam_min %.3g); point certified" % lam)
            return rec
        rec.update(verdict="PASS_REF" if better else "PASS", detail=rec.get("pass_ref_why", ""))
        return rec

    # ---- non-optimal verdicts
    if ost == "infeasible":
        ph = qpcheck.farkas_phase1(case)
        rec["phase1"] = ph
        if truth == "infeasible":
            if conv == "border" and equality_model is None:
                rec.update(verdict="BORDER", detail="")
            else:
                rec.update(verdict="PASS", detail="independent phase-1 LP violation %.3g" % ph)
        elif ph > 1e-7:
            rec.update(verdict="PASS_REF", detail="reference %s but independent phase-1 LP proves infeasible (%.3g)" % (ref["status"], ph))
        else:
            rec.update(verdict="FAIL_STATUS", detail="infeasible but truth/reference say %s; phase-1 violation %.3g" % (truth, ph))
        return rec
    if ost == "unbounded":
        if truth == "unbounded":
            rec.update(verdict="PASS", detail="")
        else:
            rec.update(verdict="FAIL_STATUS", detail="unbounded but truth/reference say %s" % truth)
        return rec
    if ost == "dual_infeasible":
        if truth in ("unbounded", "infeasible"):
            rec.update(verdict="PARTIAL", detail="dual_infeasible (recession ray verified, primal feasibility not established); truth %s" % truth)
        else:
            rec.update(verdict="FAIL_STATUS", detail="dual_infeasible but truth/reference say %s" % truth)
        return rec
    rec.update(verdict="NOANS", detail="%s: %s" % (ost, rec["msg"]))
    return rec


def work(engine, outdir, cat, i, tl):
    case = qpcase.gen(cat, i)
    d = os.path.join(outdir, "cases", cat)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, case.id + ".mps")
    with open(path, "w") as f:
        f.write(qpcase.write_mps(case))
    sol, js = path[:-4] + ".sol", path[:-4] + ".json"
    ours = run_engine(engine, path, sol, js, tl)
    x = read_sol(sol, case.n) if ours.get("status") == "optimal" else None
    rec = dict(id=case.id, cat=cat, seed=case.seed, n=case.n, m=case.m, maximize=case.maximize, qform=case.qform,
               expect=case.expect_status, note=case.note)
    try:
        rec.update(judge(case, ours, x, tl))
    except Exception:
        import traceback
        rec.update(verdict="HARNESS_ERR", detail=traceback.format_exc()[-500:])
    if rec["verdict"] in ("PASS", "PASS_REF", "BORDER", "PARTIAL"):
        for p in (path, sol, js):
            if os.path.exists(p):
                os.remove(p)
    return rec


def worker_main():
    engine, outdir, tl = sys.argv[2], sys.argv[3], float(sys.argv[4])
    proto = os.fdopen(os.dup(1), "w")  # HiGHS prints diagnostics on fd 1 (also from forked children): keep the
    os.dup2(2, 1)                      # JSON protocol on a private copy and send everything else to stderr
    for line in sys.stdin:
        cat, i = line.split()
        rec = work(engine, outdir, cat, int(i), tl)
        proto.write(json.dumps(rec, sort_keys=True, default=lambda o: None if o != o else float(o)) + "\n")
        proto.flush()


def supervise(jobs, nworkers, engine, outdir, tl, results):
    import queue
    q = queue.Queue()
    for j in jobs:
        q.put(j)
    lock = threading.Lock()
    done = [0]
    t0 = time.time()
    cap = 3 * tl + 180

    def spawn():
        return subprocess.Popen([sys.executable, os.path.abspath(__file__), "--worker", engine, outdir, str(tl)],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1)

    def loop():
        w = spawn()
        while True:
            try:
                cat, i = q.get_nowait()
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
                seed = qpcase.BASE[cat] + i
                rec = dict(id="%s-%d" % (cat, seed), cat=cat, seed=seed, verdict="HARNESS_ERR",
                           detail="worker died or hung on this case")
                try:
                    w.kill()
                except Exception:
                    pass
                w = spawn()
            with lock:
                results.append(rec)
                done[0] += 1
                if done[0] % 200 == 0:
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


def selftest(cats, limit):
    """Generator / writer / certificate self-check, no engine: the constructed optima certify, rays verify, and
    HiGHS reading the written MPS agrees with HiGHS built from the data (objective on the same point and optimum)."""
    import highspy
    bad = 0
    for cat in cats:
        n = min(limit or 10, qpcase.COUNT[cat])
        for i in range(0, n):
            case = qpcase.gen(cat, i)
            conv, lam = convexity(case)
            msgs = []
            if case.ray is not None:
                e = qpcheck.verify_ray(case)
                if e:
                    msgs.append("ray: " + e)
            if case.x0 is not None and case.expect_status == "optimal" and conv != "border" and cat not in ("psd_free_null", "ill_cond"):
                pc, kc = certify_point(case, case.x0)
                if max(pc["row_rel"], pc["bnd_rel"]) > 1e-9 or not (kc["stat_rel"] <= 1e-9):
                    msgs.append("x0 does not certify: viol %.2e %.2e stat %.2e" % (pc["row_rel"], pc["bnd_rel"], kc["stat_rel"]))
            if conv == "convex" and case.expect_status not in ("infeasible", "unbounded"):
                path = os.path.join("/tmp", "qpself_%s.mps" % case.id)
                open(path, "w").write(qpcase.write_mps(case))
                h = highspy.Highs(); h.setOptionValue("output_flag", False)
                h.readModel(path); h.run()
                st = h.modelStatusToString(h.getModelStatus())
                ref = qpcheck.reference(case, True, 30)
                if st == "Optimal" and ref["status"] == "optimal":
                    xf = np.array(h.getSolution().col_value)
                    fo = case.obj(xf)
                    o2 = ref["objective"]
                    if abs(fo - o2) > 1e-6 * max(1, abs(o2)):
                        msgs.append("file reading disagrees with data: %.10g vs %.10g (qform %s)" % (fo, o2, case.qform))
                    sg = -1 if case.maximize else 1
                    rep = sg * h.getInfo().objective_function_value  # HiGHS' own objective incl. offset, in min sense
                    if abs(rep - sg * fo) > 1e-5 * max(1, abs(fo)):
                        msgs.append("HiGHS file objective %.10g differs from data objective %.10g" % (rep, fo))
                elif st != ref["status"].capitalize() and not (st == "Optimal" and ref["status"] == "optimal"):
                    msgs.append("status via file %s vs via data %s" % (st, ref["status"]))
                os.remove(path)
            if msgs:
                bad += 1
                print("SELFTEST %s: %s" % (case.id, "; ".join(msgs)))
    print("selftest: %d problem cases" % bad)
    return bad


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        return worker_main()
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine")
    ap.add_argument("--out", default="out/adversarial_qp")
    ap.add_argument("--cats", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--time-limit", type=float, default=30.0)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    cats = [c for c in a.cats.split(",") if c] or [c for c, _, _ in qpcase.ALL]
    if a.selftest:
        sys.exit(1 if selftest(cats, a.limit) else 0)
    engine = os.path.abspath(a.engine)
    os.makedirs(a.out, exist_ok=True)
    jobs = [(c, i) for c in cats for i in range(min(a.limit, qpcase.COUNT[c]) if a.limit else qpcase.COUNT[c])]
    results = []
    t0 = time.time()
    supervise(jobs, a.jobs, engine, os.path.abspath(a.out), a.time_limit, results)
    results.sort(key=lambda r: (cats.index(r["cat"]), r["seed"]))
    with open(os.path.join(a.out, "results.jsonl"), "w") as f:
        for r in results:
            f.write(json.dumps(r, sort_keys=True, default=lambda o: None if o != o else float(o)) + "\n")
    print("ran %d cases in %.0f s -> %s" % (len(results), time.time() - t0, os.path.join(a.out, "results.jsonl")))


if __name__ == "__main__":
    main()
