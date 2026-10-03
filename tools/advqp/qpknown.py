#!/usr/bin/env python3
"""Write tests/qp_known_failures.tsv and tests/qp_repro/ from a harness run.

Every non-passing case (wrong / noans / other, plus partial, which is listed for the record) is mapped to a failure
class; the smallest case of each class (n + m) becomes a reproducer: the MPS exactly as the generator writes it, plus
an entry in expected.json that names the generator case, so tests/qp_repro replays it through the full harness checks.
Usage: qpknown.py --results out/adversarial_qp/results.jsonl [--known tests/qp_known_failures.tsv] [--repro tests/qp_repro]
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qpcase

DESC = {
    "false_optimal_gap": "reports optimal with a certified duality gap far above its own tolerance (self-reported gap ~1e-10)",
    "numfail_huge_bounds": "numerical_failure on a convex QP with finite inactive bounds of size 1e6..1e19",
    "numfail_ill_cond": "numerical_failure on a convex QP with Q condition number 1e10..1e12",
    "numfail_other": "numerical_failure on a well-posed convex QP (see detail)",
    "refuses_indefinite_q_convex_on_feasible_set": "nonconvex: Q indefinite but positive definite on the null space of the equality rows",
    "dual_infeasible_on_unbounded": "dual_infeasible (ray verified, primal feasibility not established) for a feasible unbounded model",
    "reference_error_near_psd": "HiGHS cannot solve a near-semidefinite model; the engine's point is certified (no engine fault)",
}


def classify(r):
    v, st, cat = r["verdict"], r.get("ours_status"), r["cat"]
    if v == "FAIL_KKT":
        return "false_optimal_gap"
    if st == "numerical_failure":
        return {"bounds_huge": "numfail_huge_bounds", "ill_cond": "numfail_ill_cond"}.get(cat, "numfail_other")
    if cat == "indef_convex_feasible" and st == "nonconvex":
        return "refuses_indefinite_q_convex_on_feasible_set"
    if v == "PARTIAL":
        return "dual_infeasible_on_unbounded"
    if v == "REF_UNRESOLVED":
        return "reference_error_near_psd"
    return "UNCLASSIFIED"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--known", default="tests/qp_known_failures.tsv")
    ap.add_argument("--repro", default="tests/qp_repro")
    a = ap.parse_args()
    recs = [json.loads(l) for l in open(a.results)]
    fails = [r for r in recs if r["verdict"] not in ("PASS", "PASS_REF", "BORDER")]
    best = {}
    with open(a.known, "w") as f:
        f.write("# id\tclass\treproducer\tdetail   (regenerate with tools/advqp/qpknown.py; PARTIAL rows are correct-but-weaker verdicts)\n")
        for r in sorted(fails, key=lambda r: r["id"]):
            cl = classify(r)
            rep = "-"
            if cl != "UNCLASSIFIED" and cl != "reference_error_near_psd":
                rep = "qp_%s.mps" % cl
                if cl not in best or (r["n"] + r["m"], r["seed"]) < best[cl][0]:
                    best[cl] = ((r["n"] + r["m"], r["seed"]), r)
            f.write("%s\t%s\t%s\t%s: %s\n" % (r["id"], cl, rep, r["verdict"], (r.get("detail") or "")[:140].replace("\t", " ")))
    os.makedirs(a.repro, exist_ok=True)
    exp = {}
    for cl, (_, r) in sorted(best.items()):
        cat = r["cat"]
        case = qpcase.gen(cat, r["seed"] - qpcase.BASE[cat])
        fn = "qp_%s.mps" % cl
        open(os.path.join(a.repro, fn), "w").write(qpcase.write_mps(case))
        exp[fn] = {"case": r["id"], "category": cat, "class": cl, "known_failing": True,
                   "partial_ok": False, "observed_bad": "%s: %s" % (r["verdict"], (r.get("detail") or "")[:200]),
                   "report": DESC[cl]}
    json.dump(exp, open(os.path.join(a.repro, "expected.json"), "w"), indent=1, sort_keys=True)
    print("wrote %d known-failure rows, %d reproducers" % (len(fails), len(exp)))


if __name__ == "__main__":
    main()
