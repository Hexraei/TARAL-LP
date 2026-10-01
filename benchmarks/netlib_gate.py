#!/usr/bin/env python3
"""Reproducible Netlib benchmark for the C++ engine: one command, one ledger.

Downloads the pinned public Netlib MPS files, verifies their SHA-256, obtains reference objectives
by running HiGHS on each ORIGINAL file (reference only, never part of the solve path), runs the
engine on every case, re-checks every returned solution with the independent original-model
checker (benchmarks/orig_check.py: its own parser, rows, bounds and objective incl. the
objective-row RHS constant), and writes a per-case ledger plus an environment record.

Protocols (never merged):
  60 s headline   --time-limit 60. PILOT.WE and PILOT4 are a fixed reference exclusion under this
                  protocol (HiGHS failed on them in the 60-second reference record): never passes,
                  ceiling 91/93.
  extended cap    --time-limit 300 --extended: those two are graded against the HiGHS reference.
Denominator: the 93 cases below. TRUSS is in the corpus but outside the 93 (--with-extra).
Pass rule: engine 'optimal'; HiGHS optimal on the original file; |obj - ref| <= max(1e-6, 1e-7|ref|);
checker-recomputed objective agrees within the same tolerance; relative row violation (normalised by
1 + |row rhs|) <= 1e-6; bound violation <= 1e-6. STRICT side line: both violations <= 1e-8 and
objective relative error <= 1e-8, reported beside the headline, not a pass rule.

Usage (repository root; needs g++, Python 3, NumPy and highspy):
  python benchmarks/netlib_gate.py --out results/run_60s                      # full 93 at 60 s
  python benchmarks/netlib_gate.py --out results/run_300s --time-limit 300 --extended
  python benchmarks/netlib_gate.py --out /tmp/smoke --cases afiro,e226,sc50a  # quick smoke test
"""
import argparse, csv, hashlib, json, os, platform, subprocess, sys, time, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import orig_check as oc  # independent original-model parser/checker

CORPUS_COMMIT = "56257eea85b433ce6aa67d26156b36385318fd6f"
CORPUS_URL = "https://raw.githubusercontent.com/ozy4dm/lp-data-netlib/" + CORPUS_COMMIT + "/mps_files/{}"
MANIFEST = ROOT / "benchmarks" / "netlib_manifest.sha256"
CASES_93 = """25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2 bore3d brandy
capri cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800 finnis fit1d fit1p fit2d fit2p forplan
ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7 israel kb2 lotfi maros maros-r7 modszk1 nesm perold pilot
pilot.ja pilot.we pilot4 pilot87 pilotnov recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion
scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l ship04s ship08l ship08s ship12l ship12s
sierra stair standata standgub standmps stocfor1 stocfor2 tuff vtp.base wood1p woodw""".split()
EXTRA = ["truss"]
REFERENCE_EXCLUDED_60S = {"pilot.we", "pilot4"}


def fetch_corpus(corpus):
    corpus.mkdir(parents=True, exist_ok=True)
    expected = dict(reversed(line.split()) for line in MANIFEST.read_text().splitlines() if line.strip())
    for name, digest in expected.items():
        path = corpus / name
        if not path.exists():
            urllib.request.urlretrieve(CORPUS_URL.format(name), path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            sys.exit(f"SHA-256 mismatch for {name}: corpus is not the pinned version")
    return len(expected)


def highs_reference(path, limit):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(limit))
    h.readModel(str(path))
    t = time.time()
    h.run()
    return {"status": h.modelStatusToString(h.getModelStatus()), "objective": h.getInfo().objective_function_value,
            "wall_s": time.time() - t, "highs_version": h.version()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--engine", type=Path, help="prebuilt engine binary (default: build src/ into OUT/taral)")
    ap.add_argument("--source-commit", help="commit the prebuilt --engine was built from (recorded in the summary)")
    ap.add_argument("--corpus", type=Path, default=ROOT / "corpus")
    ap.add_argument("--time-limit", type=float, default=60)
    ap.add_argument("--extended", action="store_true", help="extended-cap protocol: grade PILOT.WE/PILOT4")
    ap.add_argument("--with-extra", action="store_true", help="also run TRUSS, outside the 93")
    ap.add_argument("--cases", help="comma-separated subset (smoke test)")
    ap.add_argument("--reference-limit", type=float, default=300, help="HiGHS time limit for the reference")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    nfiles = fetch_corpus(a.corpus)

    engine = a.engine
    if engine is None:
        engine = a.out / "taral"
        srcs = sorted(str(p) for p in (ROOT / "src").glob("*.cpp"))
        subprocess.run(["g++", "-O3", "-march=native", "-std=c++17", "-o", str(engine), *srcs], check=True)

    ref_path = a.out / "highs_reference.json"
    ref = json.loads(ref_path.read_text()) if ref_path.exists() else {}
    names = a.cases.split(",") if a.cases else CASES_93 + (EXTRA if a.with_extra else [])
    protocol = f"{a.time_limit:g}s" + (" extended-cap" if a.extended else "")
    rows = []
    for n in names:
        mps = a.corpus / f"{n}.mps"
        if n not in ref:
            ref[n] = highs_reference(mps, a.reference_limit)
            ref_path.write_text(json.dumps(ref, indent=1))
        R = ref[n]
        sol, js = a.out / f"{n}.sol", a.out / f"{n}.json"
        for f in (sol, js):
            f.unlink(missing_ok=True)
        t = time.time()
        try:
            subprocess.run([str(engine), str(mps), "--time-limit", str(a.time_limit), "--sol", str(sol), "--json", str(js)],
                           timeout=a.time_limit + 15, capture_output=True)
            res = json.loads(js.read_text()) if js.exists() else {"status": "no_output"}
        except subprocess.TimeoutExpired:
            res = {"status": "killed_after_cap"}
        row = {"case": n, "in_denominator": n in CASES_93, "engine_status": res.get("status"),
               "objective": res.get("objective"), "iterations": res.get("iterations"), "wall_s": round(time.time() - t, 3),
               "message": res.get("message", ""), "reference_status": R["status"], "reference_objective": R["objective"],
               "verdict": "", "pass": False, "strict": False}
        if n in REFERENCE_EXCLUDED_60S and not a.extended:
            row["verdict"] = "reference_excluded_60s"
        elif R["status"] != "Optimal":
            row["verdict"] = "reference_not_optimal"
        elif res.get("status") == "optimal" and sol.exists() and res.get("objective") is not None:
            x = {}
            for line in sol.read_text().splitlines():
                j, v = line.rsplit(None, 1)
                x[j] = float(v)
            chk = oc.check(oc.parse_mps(mps), x)
            ho, ro = R["objective"], res["objective"]
            tol = max(1e-6, 1e-7 * abs(ho))
            ok = (abs(ro - ho) <= tol and abs(chk["objective"] - ro) <= tol
                  and chk["row_violation_rel"] <= 1e-6 and chk["bound_violation"] <= 1e-6)
            row.update(checker_objective=chk["objective"], abs_objective_error=abs(ro - ho),
                       rel_objective_error=abs(ro - ho) / max(1.0, abs(ho)), row_violation_rel=chk["row_violation_rel"],
                       bound_violation=chk["bound_violation"], worst_row=chk["worst_row"], pass_=ok)
            row["pass"] = ok
            row["strict"] = ok and chk["row_violation_rel"] <= 1e-8 and chk["bound_violation"] <= 1e-8 \
                and abs(ro - ho) <= 1e-8 * max(1.0, abs(ho))
            row["verdict"] = "pass" if ok else "WRONG_OR_INFEASIBLE_ANSWER"
            row.pop("pass_")
        else:
            row["verdict"] = res.get("status") or "no_valid_result"
        rows.append(row)
        print(f"{n:10} {row['verdict']:28} {row['wall_s']:7.1f}s", flush=True)

    fields = ["case", "in_denominator", "verdict", "pass", "strict", "engine_status", "objective", "reference_status",
              "reference_objective", "checker_objective", "abs_objective_error", "rel_objective_error", "row_violation_rel",
              "bound_violation", "worst_row", "iterations", "wall_s", "message"]
    with open(a.out / "ledger.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    core = [r for r in rows if r["in_denominator"]]
    passes = [r for r in core if r["pass"]]
    wrong = [r["case"] for r in core if r["verdict"] == "WRONG_OR_INFEASIBLE_ANSWER"]
    if a.engine is None:  # built from the working tree here
        commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "src"], capture_output=True, text=True).stdout)
    else:
        commit, dirty = a.source_commit or "unknown (prebuilt binary)", None
    summary = {"protocol": protocol, "cases_run": len(core), "denominator": 93, "passes": len(passes),
               "strict_passes": sum(r["strict"] for r in passes), "wrong_answers": wrong,
               "failures": {r["case"]: r["verdict"] for r in core if not r["pass"]},
               "extra_outside_denominator": {r["case"]: r["verdict"] for r in rows if not r["in_denominator"]},
               "engine_binary": str(engine), "repo_commit": commit, "src_has_uncommitted_changes": dirty,
               "corpus": f"ozy4dm/lp-data-netlib@{CORPUS_COMMIT} ({nfiles} files, SHA-256 verified)",
               "tolerances": {"objective": "max(1e-6, 1e-7|ref|)", "row_rel": 1e-6, "bound": 1e-6, "strict": 1e-8},
               "environment": {"platform": platform.platform(), "processor": next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")),
                                                 platform.processor()) if os.path.exists("/proc/cpuinfo") else platform.processor(),
                               "cpu_count": os.cpu_count(), "python": platform.python_version(),
                               "compiler": subprocess.run(["g++", "--version"], capture_output=True, text=True).stdout.split("\n")[0]}}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(f"\nPROTOCOL {protocol} | PASSES {len(passes)}/93 (of {len(core)} run) | STRICT {summary['strict_passes']}/"
          f"{len(passes)} passes | WRONG {wrong} | extra {summary['extra_outside_denominator']}")


if __name__ == "__main__":
    main()
