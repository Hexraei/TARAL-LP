"""Run existing small benchmarks and record a reproducible CPU smoke result."""
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

import highspy
import orig_check
from mps_semantics_check import CASES

ROOT = Path(__file__).resolve().parents[1]
PINS = {"numpy": "2.2.6", "highspy": "1.11.0"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lp_methods(engine):
    """Check each LP method on an existing fixture, not a new model generator."""
    text, status, objective, _ = CASES["free_format_d_exponent"]
    with tempfile.TemporaryDirectory(prefix="taral_lp_smoke_") as tmp:
        mps = Path(tmp) / "lp.mps"
        mps.write_text(text)
        model = orig_check.parse_mps(mps)
        ref = highspy.Highs()
        ref.setOptionValue("output_flag", False)
        if ref.readModel(str(mps)) != highspy.HighsStatus.kOk:
            raise RuntimeError("HiGHS could not read LP fixture")
        ref.run()
        if ref.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            raise RuntimeError("HiGHS did not solve LP fixture")
        for method in ("simplex", "dual", "ipm"):
            js, sol = Path(tmp) / "result.json", Path(tmp) / "result.sol"
            js.unlink(missing_ok=True)
            sol.unlink(missing_ok=True)
            subprocess.run([str(engine), str(mps), "--method", method,
                            "--time-limit", "10", "--json", str(js), "--sol", str(sol)],
                           check=True, capture_output=True, text=True, timeout=30)
            result = json.loads(js.read_text())
            if result["status"] != status:
                raise RuntimeError(f"{method}: expected {status}, got {result}")
            point = dict((name, float(value)) for name, value in
                         (line.split() for line in sol.read_text().splitlines()))
            if set(point) != set(model["corder"]):
                raise RuntimeError(f"{method}: incomplete solution columns")
            if not all(math.isfinite(v) for v in point.values()):
                raise RuntimeError(f"{method}: non-finite solution")
            checked = orig_check.check(model, point)
            for want in (objective, ref.getInfo().objective_function_value, result["objective"]):
                if want is None or not math.isfinite(want) or abs(checked["objective"] - want) > 1e-6 * max(1, abs(want)):
                    raise RuntimeError(f"{method}: objective mismatch {checked}, {want}")
            if checked["row_violation_rel"] > 1e-6 or checked["bound_violation"] > 1e-6:
                raise RuntimeError(f"{method}: infeasible point {checked}")
            print(f"PASS LP {method}: objective {checked['objective']:.10g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=57463)
    args = ap.parse_args()
    engine, out = args.engine.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    versions = {name: importlib.metadata.version(name) for name in PINS}
    if versions != PINS:
        ap.error(f"use requirements-smoke.txt; installed {versions}, expected {PINS}")
    if not engine.is_file() or not os.access(engine, os.X_OK):
        ap.error(f"engine is not executable: {engine}")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                         capture_output=True, text=True)
    cache = engine.parent / "CMakeCache.txt"
    build_settings = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            if re.match(r"CMAKE_(CXX_COMPILER|CXX_FLAGS|CXX_FLAGS_RELEASE|BUILD_TYPE):", line):
                key, value = line.split("=", 1)
                build_settings[key.split(":", 1)[0]] = value
    compiler = build_settings.get("CMAKE_CXX_COMPILER")
    compiler_version = None
    if compiler:
        version = subprocess.run([compiler, "--version"], capture_output=True, text=True, timeout=10)
        compiler_version = version.stdout.splitlines()[0] if version.returncode == 0 else None
    metadata = {
        "schema_version": 1, "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "seed": args.seed, "per_suite": 1, "workers": 1,
        "python": sys.version, "platform": platform.platform(), "dependencies": versions,
        "build_settings": build_settings, "compiler_version": compiler_version,
        "engine": str(engine), "engine_sha256": sha256(engine),
        "source_commit": git.stdout.strip() if git.returncode == 0 else None,
        "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sorted((ROOT / "src").glob("*")) if p.is_file()},
        "refinery_sha256": sha256(ROOT / "examples/refinery_williams.py"),
        "requirements_sha256": sha256(ROOT / "benchmarks/requirements-smoke.txt"),
        "benchmark_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in
                             sorted((ROOT / "benchmarks").glob("*.py"))},
        "checks": [],
    }
    start = time.monotonic()
    try:
        lp_methods(engine)
        metadata["checks"].append({"name": "lp_methods", "passed": True})
    except Exception as exc:
        metadata["checks"].append({"name": "lp_methods", "passed": False, "error": str(exc)})
        print(f"FAIL LP methods: {exc}")
    suites = [
        ("mps", "benchmarks/mps_semantics_check.py", []),
        ("milp", "benchmarks/milp_tests.py", ["--seed", str(args.seed), "--per-suite", "1", "--workers", "1"]),
        ("qp", "benchmarks/qp_tests.py", ["--seed", str(args.seed), "--per-suite", "1", "--workers", "1"]),
        ("refinery", "examples/refinery_williams.py", []),
    ]
    for name, script, extra in suites:
        command = [sys.executable, str(ROOT / script), "--engine", str(engine), *extra]
        entry = {"name": name, "command": command, "passed": False, "log": str(out / f"{name}.log")}
        try:
            run = subprocess.run(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, timeout=180)
            (out / f"{name}.log").write_text(run.stdout)
            entry.update(returncode=run.returncode, passed=run.returncode == 0)
            # Existing MILP/QP benchmarks only fail on wrong, not undecided.
            # Keep the deliberately node-limited MILP case, but reject any other undecided case.
            if name in ("milp", "qp"):
                totals = re.search(r"TOTAL (\d+) cases: ok (\d+), undecided (\d+), wrong (\d+)", run.stdout)
                if totals is None:
                    entry["passed"] = False
                    entry["error"] = "missing benchmark summary"
                else:
                    entry["tally"] = dict(zip(("total", "ok", "undecided", "wrong"), map(int, totals.groups())))
                    if entry["tally"]["total"] != 8 or entry["tally"]["wrong"]:
                        entry["passed"] = False
                        entry["error"] = "unexpected case count or wrong result"
                    unknown = [line for line in run.stdout.splitlines() if line.startswith("UNDECIDED")]
                    expected = name == "milp" and len(unknown) == entry["tally"]["undecided"] and all(re.match(r"UNDECIDED\s+h_limits\s+", line) for line in unknown)
                    if entry["tally"]["undecided"] and not expected:
                        entry["passed"] = False
                        entry["error"] = "unexpected undecided case"
            print(f"{'PASS' if entry['passed'] else 'FAIL'} {name}: {entry['log']}")
        except Exception as exc:
            entry["error"] = str(exc)
            print(f"FAIL {name}: {exc}")
        metadata["checks"].append(entry)
    metadata["wall_s"] = time.monotonic() - start
    metadata["passed"] = all(check["passed"] for check in metadata["checks"])
    (out / "summary.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"{'PASS' if metadata['passed'] else 'FAIL'} smoke: {out / 'summary.json'}")
    return 0 if metadata["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
