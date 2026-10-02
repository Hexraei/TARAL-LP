#!/usr/bin/env python3
"""Pinned small-MIPLIB ledger for the C++ MILP solver (local measurement, not a published result).

Source: the official MIPLIB 3 archive (ZIB), pinned by sha256; every extracted instance file is
hashed too. Reference optima: the published MIPLIB 3 catalogue value ("INT SOLN", all marked
optimal for the instances used here) and a HiGHS solve of the same file (cached in
out/miplib/highs_ref.json). Every taral solution is re-checked against the ORIGINAL model as read
by HiGHS (rows, bounds, integrality, objective), so the check shares no code with taral's reader.

Pass rules (minimisation sense; sign flipped for MAX models):
  solved   status optimal, point feasible, |obj - ref| <= 1e-6 * max(1, |ref|)
  limit    time_limit/node_limit and honest: bound <= ref (+tol) and incumbent (if any) >= ref (-tol)
  FAIL     anything else (wrong optimum, infeasible point, invalid bound, wrong status, crash)

Usage: ~/.venvs/taral-gpu/bin/python benchmarks/milp_miplib.py --time-limit 60 [--engine out/taral]
Writes out/miplib/ledger_<T>s.csv and out/miplib/sources.csv.
"""
import argparse, csv, hashlib, json, math, os, subprocess, tarfile, time, urllib.request

import numpy as np
import highspy

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DIR = os.path.join(ROOT, "out", "miplib")
URL = "https://miplib2010.zib.de/miplib3/miplib3.tar.gz"
SHA256 = "f07be4c2b3cad25b23b0c3085ec98412524989dbd24f4d385e5ffa1e3e4c118f"
# instance -> published MIPLIB 3 optimum (miplib.cat, "INT SOLN")
PUBLISHED = {
    "p0033": 3089, "p0201": 7615, "p0282": 258411, "p0548": 8691, "egout": 568.101, "flugpl": 1201500,
    "lseu": 1120, "mod008": 307, "stein27": 18, "stein45": 30, "bell5": 8966406.49, "bell3a": 878430.32,
    "misc03": 3360, "misc07": 2810, "enigma": 0.0, "gt2": 21166.0, "khb05250": 106940226, "vpm1": 20,
    "vpm2": 13.75, "dcmulti": 188182, "rgn": 82.1999, "pk1": 11.0, "pp08a": 7350.0, "noswot": -43,
    "fixnet6": 3983, "blend2": 7.598985,
}
# published values above are rounded in the catalogue; compare with this relative tolerance
PUBLISHED_TOL = 1e-6


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def fetch():
    os.makedirs(DIR, exist_ok=True)
    tar = os.path.join(DIR, "miplib3.tar.gz")
    if not os.path.exists(tar):
        urllib.request.urlretrieve(URL, tar)
    if sha(tar) != SHA256:
        raise SystemExit("archive sha256 mismatch: %s" % sha(tar))
    with tarfile.open(tar) as t:  # only the pinned instances
        t.extractall(DIR, members=[t.getmember("miplib3/" + n) for n in PUBLISHED], filter="data")
    for n in PUBLISHED:  # HiGHS picks the reader by extension
        os.replace(os.path.join(DIR, "miplib3", n), os.path.join(DIR, "miplib3", n + ".mps"))
    with open(os.path.join(DIR, "sources.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "source", "sha256"])
        w.writerow(["miplib3.tar.gz", URL, SHA256])
        for n in PUBLISHED:
            p = os.path.join(DIR, "miplib3", n + ".mps")
            w.writerow(["miplib3/%s.mps" % n, URL + " (member miplib3/%s)" % n, sha(p)])


def highs_model(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(path)
    return h


def highs_ref(names, limit):
    cache_path = os.path.join(DIR, "highs_ref.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    for n in names:
        if n in cache:
            continue
        h = highs_model(os.path.join(DIR, "miplib3", n + ".mps"))
        h.setOptionValue("mip_rel_gap", 1e-9)
        h.setOptionValue("mip_abs_gap", 1e-9)
        h.setOptionValue("time_limit", float(limit))
        t0 = time.time()
        h.run()
        info = h.getInfo()
        cache[n] = {"status": h.modelStatusToString(h.getModelStatus()), "objective": info.objective_function_value,
                    "bound": info.mip_dual_bound, "wall_s": time.time() - t0, "highs_version": h.version()}
        json.dump(cache, open(cache_path, "w"), indent=1)
    return cache


def check_point(path, sol):
    """Max violation (rows/bounds relative to 1+|b|, integrality absolute) and objective, model read by HiGHS."""
    h = highs_model(path)
    lp = h.getLp()
    vals = {}
    for line in open(sol):
        name, v = line.rsplit(None, 1)
        vals[name.strip()] = float(v)
    x = np.array([vals[c] for c in lp.col_names_])
    a = lp.a_matrix_
    start, index, value = np.array(a.start_), np.array(a.index_), np.array(a.value_)
    act = np.zeros(lp.num_row_)
    for j in range(lp.num_col_):
        s, e = start[j], start[j + 1]
        act[index[s:e]] += value[s:e] * x[j]
    worst = 0.0
    for v, lo, up in list(zip(act, lp.row_lower_, lp.row_upper_)) + list(zip(x, lp.col_lower_, lp.col_upper_)):
        if lo > -1e30:
            worst = max(worst, (lo - v) / (1 + abs(lo)))
        if up < 1e30:
            worst = max(worst, (v - up) / (1 + abs(up)))
    integ = [int(t) for t in lp.integrality_] if len(lp.integrality_) else [0] * lp.num_col_
    for j, t in enumerate(integ):
        if t:
            worst = max(worst, abs(x[j] - round(x[j])))
    obj = lp.offset_ + float(np.dot(lp.col_cost_, x))
    sense = -1 if str(lp.sense_).endswith("kMaximize") else 1
    return worst, obj, sense


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", default=os.path.join(ROOT, "out", "taral"))
    ap.add_argument("--time-limit", type=float, default=60)
    ap.add_argument("--highs-limit", type=float, default=600)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    fetch()
    names = a.only or list(PUBLISHED)
    ref = highs_ref(names, a.highs_limit)
    out = os.path.join(DIR, "ledger_%gs.csv" % a.time_limit)
    cols = ["instance", "rows", "cols", "ints", "status", "objective", "best_bound", "gap", "nodes", "wall_s",
            "ref_published", "ref_highs", "highs_status", "abs_err", "rel_err", "max_violation", "verdict", "reason"]
    rows = []
    for n in names:
        path = os.path.join(DIR, "miplib3", n + ".mps")
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
            r = {"status": "crash", "objective": None, "best_bound": None, "gap": None, "nodes": None,
                 "wall_s": time.time() - t0, "message": repr(e)}
        lp = highs_model(path).getLp()
        ints = sum(1 for t in lp.integrality_ if int(t)) if len(lp.integrality_) else 0
        pub, hr = PUBLISHED[n], ref[n]
        refv = hr["objective"] if hr["status"] == "Optimal" else pub
        viol, abs_err, rel_err, verdict, reason = None, None, None, "FAIL", ""
        sense = 1
        if os.path.exists(sol):
            viol, obj, sense = check_point(path, sol)
            if r["objective"] is not None and abs(obj - r["objective"]) > 1e-6 * max(1, abs(obj)):
                reason = "reported objective %r != recomputed %r; " % (r["objective"], obj)
        tol = 1e-6 * max(1, abs(refv))
        if hr["status"] == "Optimal" and abs(hr["objective"] - pub) > PUBLISHED_TOL * max(1, abs(pub)) + 0.01:
            reason += "HiGHS %r differs from published %r; " % (hr["objective"], pub)
        if r["objective"] is not None:
            abs_err = abs(r["objective"] - refv)
            rel_err = abs_err / max(1, abs(refv))
        st = r["status"]
        if viol is not None and viol > 1e-6:
            reason += "point violates the original model by %.3g" % viol
        elif st == "optimal":
            if abs_err is not None and abs_err <= tol:
                verdict = "solved"
            else:
                reason += "optimal claimed with objective %r vs reference %r" % (r["objective"], refv)
        elif st in ("time_limit", "node_limit"):
            b = r["best_bound"]
            if b is not None and sense * b > sense * refv + tol:
                reason += "bound %r cuts off the reference %r" % (b, refv)
            elif r["objective"] is not None and sense * r["objective"] < sense * refv - tol:
                reason += "incumbent %r better than the reference %r" % (r["objective"], refv)
            else:
                verdict = "limit"
                reason += "not solved within %gs (%s)" % (a.time_limit, "incumbent" if r["objective"] is not None
                                                          else "no incumbent")
        else:
            reason += "status %s: %s" % (st, r.get("message", ""))
        rows.append([n, lp.num_row_, lp.num_col_, ints, st, r["objective"], r["best_bound"], r["gap"], r["nodes"],
                     round(r["wall_s"], 3), pub, hr["objective"], hr["status"], abs_err, rel_err, viol, verdict,
                     reason.strip("; ")])
        print(" ".join(str(v) for v in rows[-1]), flush=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
    tally = {}
    for row in rows:
        tally.setdefault(row[-2], []).append(row[0])
    for k, v in sorted(tally.items()):
        print("%s %d: %s" % (k, len(v), " ".join(v)))
    print("ledger:", out)


if __name__ == "__main__":
    main()
