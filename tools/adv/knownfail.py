#!/usr/bin/env python3
"""Turn failing records of a harness run into lines of tests/known_failures.tsv.

Each failing case is mapped to a failure class and, where the engine is the one that is wrong, to the reproducer in
tests/repro/ that isolates the class. Cases that no rule covers are written as UNCLASSIFIED so they show up in review
rather than being silently accepted.

Usage: knownfail.py --results out/adversarial/results.jsonl [--cats a,b] >> tests/known_failures.tsv
Columns: id <TAB> class <TAB> reproducer-or-dash <TAB> verdict and one-line detail
"""
import argparse
import json


def classify(r):
    cat, v, o, ref = r["cat"], r["verdict"], r.get("ours_status"), r.get("ref_status")
    det = r.get("detail") or ""
    if cat == "coef_range":
        if o == "infeasible":
            return "coef_range_false_infeasible", "coef_range_false_infeasible.mps"
        if o == "unbounded":
            return "coef_range_false_unbounded", "coef_range_false_unbounded.mps"
        if o in ("time_limit", "iteration_limit"):
            return "coef_range_stall", "coef_range_stall.mps"
        if o == "numerical_failure":
            return "coef_range_numfail", "coef_range_numfail.mps"
    if cat == "huge_bounds":
        if o == "optimal" and ref == "unbounded":
            return "huge_bounds_sentinel_unbounded", "huge_bounds_sentinel_unbounded.mps"
        if o == "numerical_failure":
            return "huge_bounds_sentinel_numfail", "huge_bounds_sentinel_numfail.mps"
    if cat == "near_singular":
        if o == "infeasible":
            return "near_singular_false_infeasible", "near_singular_false_infeasible.mps"
        if o == "optimal":
            return "near_singular_objective", "near_singular_objective.mps"
        if o in ("time_limit", "iteration_limit"):
            return "near_singular_stall", "near_singular_stall.mps"
    if cat == "big_values" and o in ("numerical_failure", "time_limit"):
        return "big_values_numfail", "big_values_numfail.mps"
    if v == "REF_UNRESOLVED":
        if ref in ("ref_hang", "ref_crash") or (ref or "").startswith("other:"):
            return "reference_solver_defect", "-"
        if ref == "time_limit" and o == "time_limit":
            return "both_time_limit", "-"
    if v == "HARNESS_ERR":
        return "harness_worker_died", "-"
    return "UNCLASSIFIED", "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--cats", default="")
    a = ap.parse_args()
    cats = set(c for c in a.cats.split(",") if c)
    for ln in open(a.results):
        r = json.loads(ln)
        if r["verdict"] in ("PASS", "PASS_REF") or (cats and r["cat"] not in cats):
            continue
        k, rep = classify(r)
        print("%s\t%s\t%s\t%s: %s" % (r["id"], k, rep, r["verdict"], (r.get("detail") or "").replace("\t", " ")[:140]))


if __name__ == "__main__":
    main()
