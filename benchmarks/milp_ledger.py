#!/usr/bin/env python3
"""Before/after ledger for MILP engine changes (measurement tool; HiGHS is the reference only).

  random  the cpp-engine/tools/milp_check.py generator, N cases per seed, solved by --engine at --time-limit
          and graded against HiGHS (scipy.optimize.milp on the generator data) plus an independent
          feasibility/objective check of the engine's point. Reference results are cached per (seed, case).
  miplib  the pinned small-MIPLIB list of benchmarks/milp_miplib.py (same checker, same reference cache).
  compare two ledgers of the same kind joined per instance; losses are listed first.

Every ledger is a CSV plus a .meta.json that records the exact command lines, seeds, commit and tool versions.
Cloud timings are indicative only; node counts are deterministic unless a run hits the time limit.

  python benchmarks/milp_ledger.py random --engine OUT/taral --commit HASH --out DIR/random_before.csv
  python benchmarks/milp_ledger.py miplib --engine OUT/taral --commit HASH --out DIR/miplib_before.csv
  python benchmarks/milp_ledger.py compare DIR/random_before.csv DIR/random_after.csv
"""
import argparse, csv, json, os, platform, random, subprocess, sys, tempfile, time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SEEDS = [26119, 3, 7, 26120, 31, 47]
CHECK = os.path.join(ROOT, "cpp-engine", "tools", "milp_check.py")
FIELDS = ["family", "seed", "instance", "rows", "cols", "ints", "status", "objective", "best_bound", "nodes",
          "lp_iterations", "wall_s", "prop_tightened", "prop_crossed", "prop_crossed_lp_infeasible", "prop_pruned", "rc_fixed", "rc_skipped", "audit", "ref_status",
          "ref_objective", "verdict", "detail"]
COUNTERS = ["prop_tightened", "prop_crossed", "prop_crossed_lp_infeasible", "prop_pruned", "rc_fixed", "rc_skipped"]  # absent in engines that predate them


def load_check():
    """milp_check.py runs main() on import; load its helpers without running it."""
    src = open(CHECK).read().rsplit("\nmain()", 1)[0]
    ns = {"__name__": "milp_check_helpers"}
    exec(compile(src, CHECK, "exec"), ns)
    return ns


def meta(a, cmds, extra=None):
    import highspy, scipy
    m = {"commit": a.commit, "engine": os.path.abspath(a.engine), "time_limit_s": a.time_limit,
         "command_lines": cmds, "tool_command": " ".join(sys.argv), "date": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
         "platform": platform.platform(), "processor": platform.processor(), "python": sys.version.split()[0],
         "scipy": scipy.__version__, "highspy": getattr(highspy, "__version__", "?"), "cpus": os.cpu_count(),
         "note": "cloud container; timings are indicative only"}
    m.update(extra or {})
    json.dump(m, open(a.out.rsplit(".", 1)[0] + ".meta.json", "w"), indent=1)


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})


def run_random(a):
    ck = load_check()
    cache_path = a.ref_cache
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    work = tempfile.mkdtemp(prefix="milp_ledger_")
    rows, cmds = [], []
    for seed in a.seeds:
        rng = random.Random(seed)
        for k in range(a.n):
            g = ck["gen"](rng, k)
            key = "%d/%d" % (seed, k)
            p = os.path.join(work, "s%d_%s.mps" % (seed, g["name"]))
            ck["write_mps"](p, **g)
            if key not in cache:
                st, ob = ck["ref_direct"](g)
                cache[key] = [st, ob]
            hs, ho = cache[key]
            js, sp = p + ".json", p + ".sol"
            cmd = [a.engine, p, "--time-limit", str(a.time_limit), "--json", js, "--sol", sp]
            if k == 0:
                cmds.append("seed %d: %s  (model written by milp_check.write_mps from gen(random.Random(%d), k), k=0..%d)"
                            % (seed, " ".join(["taral", "CASE.mps", "--time-limit", str(a.time_limit), "--json", "CASE.json",
                                               "--sol", "CASE.sol"]), seed, a.n - 1))
            t0 = time.time()
            try:
                subprocess.run(cmd, check=False, capture_output=True, timeout=a.time_limit * 3 + 30)
                r = json.load(open(js))
            except Exception as e:
                r = {"status": "crash", "message": repr(e), "objective": None, "best_bound": None, "nodes": None,
                     "iterations": None, "wall_s": time.time() - t0}
            ref_cls = "optimal" if hs == "Optimal" else ("infeasible" if hs == "Infeasible" else hs)
            verdict, detail = "match", ""
            if r["status"] == ref_cls:
                if ref_cls == "optimal":
                    feas, o = ck["verify_sol"](g, sp) if os.path.exists(sp) else (False, None)
                    if not feas or abs(o - r["objective"]) > 1e-6 * max(1, abs(o)):
                        verdict, detail = "WRONG", "reported optimum is not a verified feasible point with that objective"
                    elif abs(r["objective"] - ho) > 1e-6 * max(1, abs(ho)):
                        better = r["objective"] > ho if g["maxi"] else r["objective"] < ho
                        verdict, detail = ("ref_side_suboptimal", "verified taral point better than reference %r" % ho) \
                            if better else ("WRONG", "objective %r vs reference %r" % (r["objective"], ho))
            elif r["status"] == "optimal" and ref_cls == "infeasible":
                feas, o = ck["verify_sol"](g, sp) if os.path.exists(sp) else (False, None)
                verdict, detail = ("ref_side_infeasible", "reference infeasible; taral point verified feasible") \
                    if feas and abs(o - r["objective"]) <= 1e-6 * max(1, abs(o)) else ("WRONG", "optimal claimed, point fails")
            elif r["status"] in ("time_limit", "node_limit"):
                verdict, detail = "unsolved", r["status"]
            else:
                verdict, detail = "WRONG", "taral %s vs reference %s %r %s" % (r["status"], hs, ho, r.get("message", ""))
            rows.append(dict(family="random", seed=seed, instance=g["name"], rows=len(g["rows"]), cols=len(g["cols"]),
                             ints=sum(1 for x in g["ints"] if x), status=r["status"], objective=r.get("objective"),
                             best_bound=r.get("best_bound"), nodes=r.get("nodes"), lp_iterations=r.get("iterations"),
                             wall_s=round(r["wall_s"], 4), ref_status=hs, ref_objective=ho, verdict=verdict, detail=detail,
                             audit=json.dumps(r["audit"]) if "audit" in r else "", **{c: r.get(c, "") for c in COUNTERS}))
        json.dump(cache, open(cache_path, "w"))
    write_csv(a.out, rows)
    meta(a, cmds, {"seeds": a.seeds, "cases_per_seed": a.n, "reference": "scipy.optimize.milp (HiGHS) on generator data"})
    summarize(rows)


def run_miplib(a):
    sys.path.insert(0, os.path.join(ROOT, "benchmarks"))
    import milp_miplib as mm
    cmd = [sys.executable, os.path.join(ROOT, "benchmarks", "milp_miplib.py"), "--engine", a.engine,
           "--time-limit", str(a.time_limit)] + (["--only"] + a.only if a.only else [])
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL if a.quiet else None)
    src = os.path.join(mm.DIR, "ledger_%gs.csv" % a.time_limit)
    rows = []
    for r in csv.DictReader(open(src)):
        # wrong = a verdict of FAIL (the checker's own rule); limit = honest time/node limit; solved
        verdict = {"solved": "match", "limit": "unsolved"}.get(r["verdict"], "WRONG")
        try:
            ej = json.load(open(os.path.join(mm.DIR, r["instance"] + ".json")))
        except Exception:
            ej = {}
        rows.append(dict(family="miplib", seed="", instance=r["instance"], rows=r["rows"], cols=r["cols"], ints=r["ints"],
                         lp_iterations=ej.get("iterations", ""), **{c: ej.get(c, "") for c in COUNTERS},
                         status=r["status"], objective=r["objective"], best_bound=r["best_bound"], nodes=r["nodes"],
                         wall_s=r["wall_s"], ref_status=r["highs_status"], ref_objective=r["ref_highs"],
                         verdict=verdict, detail=r["reason"]))
    write_csv(a.out, rows)
    meta(a, [" ".join(cmd)], {"instances": [r["instance"] for r in rows], "miplib_archive_sha256": mm.SHA256,
                              "published_optima": mm.PUBLISHED})
    summarize(rows)


def summarize(rows):
    tally = {}
    for r in rows:
        tally[r["verdict"]] = tally.get(r["verdict"], 0) + 1
    solved = sum(1 for r in rows if r["status"] in ("optimal", "infeasible"))
    nodes = sum(int(r["nodes"] or 0) for r in rows)
    cnt = {c: sum(int(r.get(c) or 0) for r in rows) for c in COUNTERS}
    if any(r.get("audit") for r in rows):  # [nodes, LP infeasible, LP feasible, LP other, int-empty, int-feasible, undecided, cutoff-only]
        tot = [0] * 10
        for r in rows:
            for k, v in enumerate(json.loads(r["audit"]) if r.get("audit") else []): tot[k] += v
        cnt["audit[nodes,lp_infeasible,lp_feasible,lp_other,int_empty,int_feasible,undecided,cutoff_only,rc_checked,rc_bad]"] = tot
    print("cases=%d solved(status optimal/infeasible)=%d total_nodes=%d verdicts=%s counters=%s" % (len(rows), solved, nodes, tally, cnt))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def compare(a):
    B = {(r["seed"], r["instance"]): r for r in csv.DictReader(open(a.before))}
    A = {(r["seed"], r["instance"]): r for r in csv.DictReader(open(a.after))}
    assert B.keys() == A.keys(), "ledgers cover different instances"
    losses, gains, same_nodes = [], [], 0
    tot = {"b_nodes": 0, "a_nodes": 0, "b_t": 0.0, "a_t": 0.0, "b_solved": 0, "a_solved": 0}
    both_solved = {"b_nodes": 0, "a_nodes": 0, "b_t": 0.0, "a_t": 0.0, "n": 0}
    lines = []
    for k in sorted(B, key=lambda k: (str(k[0]), k[1])):
        b, r = B[k], A[k]
        bs, as_ = b["status"] in ("optimal", "infeasible"), r["status"] in ("optimal", "infeasible")
        tot["b_solved"] += bs; tot["a_solved"] += as_
        for s, x in (("b", b), ("a", r)):
            tot[s + "_nodes"] += int(x["nodes"] or 0); tot[s + "_t"] += float(x["wall_s"] or 0)
        if bs and as_:
            for s, x in (("b", b), ("a", r)):
                both_solved[s + "_nodes"] += int(x["nodes"] or 0); both_solved[s + "_t"] += float(x["wall_s"] or 0)
            both_solved["n"] += 1
        bo, ao = num(b["objective"]), num(r["objective"])
        obj_changed = (bo is None) != (ao is None) or (bo is not None and abs(bo - ao) > 1e-6 * max(1, abs(bo)))
        if b["status"] != r["status"] or obj_changed or b["verdict"] != r["verdict"]:
            tag = "LOSS" if (r["verdict"] == "WRONG" and b["verdict"] != "WRONG") or (bs and not as_) else \
                  "GAIN" if (as_ and not bs) or (b["verdict"] == "WRONG" and r["verdict"] != "WRONG") else "CHANGE"
            (losses if tag == "LOSS" else gains).append(k)
            lines.append("%-6s %-6s %-9s status %s -> %s  obj %s -> %s  verdict %s -> %s  nodes %s -> %s  t %s -> %s"
                         % (tag, k[0], k[1], b["status"], r["status"], b["objective"], r["objective"], b["verdict"],
                            r["verdict"], b["nodes"], r["nodes"], b["wall_s"], r["wall_s"]))
        same_nodes += b["nodes"] == r["nodes"]
    print("instances:", len(B), " unchanged node counts:", same_nodes)
    print("solved (optimal/infeasible):  before %d  after %d" % (tot["b_solved"], tot["a_solved"]))
    print("total nodes (all):            before %d  after %d" % (tot["b_nodes"], tot["a_nodes"]))
    print("total wall s (all):           before %.1f  after %.1f  (indicative)" % (tot["b_t"], tot["a_t"]))
    print("instances solved by both: %d  nodes before %d after %d  wall before %.1f after %.1f" % (
        both_solved["n"], both_solved["b_nodes"], both_solved["a_nodes"], both_solved["b_t"], both_solved["a_t"]))
    print("verdict WRONG: before %d  after %d" % (sum(x["verdict"] == "WRONG" for x in B.values()),
                                                  sum(x["verdict"] == "WRONG" for x in A.values())))
    print("losses: %d  gains/changes: %d" % (len(losses), len(gains)))
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("random", "miplib"):
        p = sub.add_parser(name)
        p.add_argument("--engine", required=True)
        p.add_argument("--commit", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--time-limit", type=float, default=60)
        if name == "random":
            p.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
            p.add_argument("--n", type=int, default=100)
            p.add_argument("--ref-cache", default=os.path.join(ROOT, "out", "milp_random_ref.json"))
        else:
            p.add_argument("--only", nargs="*")
            p.add_argument("--quiet", action="store_true")
    c = sub.add_parser("compare")
    c.add_argument("before")
    c.add_argument("after")
    a = ap.parse_args()
    os.makedirs(os.path.join(ROOT, "out"), exist_ok=True)
    {"random": run_random, "miplib": run_miplib, "compare": compare}[a.mode](a)


if __name__ == "__main__":
    main()
