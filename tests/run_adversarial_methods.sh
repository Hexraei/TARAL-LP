#!/usr/bin/env bash
# Adversarial suite across solver paths: the default primal-simplex path ("before") and the dual simplex (--method dual)
# and interior-point (--method ipm) paths ("after"), on the base suite plus the big-model categories (tools/adv/cats/big_*.py,
# 500-3000 rows, seeds 59000-64011). One engine build, same cases, same gates and the same 60 s per-case cap for every path.
#   tests/run_adversarial_methods.sh                 everything (about 1.5 hours on 4 cores)
#   QUICK=1 tests/run_adversarial_methods.sh         10 cases per category
#   TARAL=/path/to/engine JOBS=8 OUT=dir TL=60 tests/run_adversarial_methods.sh
# The --method flag is only honoured for pure LPs (src/main.cpp routes every MILP to branch and bound), so dual and ipm run
# the LP categories; simplex runs all categories. Every stage's exit status is recorded and printed, none is masked.
# Exit status: 0 when no path has a NEW failure (a case that passes on simplex and fails on dual/ipm) and the simplex path has
# only known failures, 1 otherwise (including a failing reproducer replay), 2 when the suite could not run.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$ROOT/out/methods}
JOBS=${JOBS:-4}
TL=${TL:-60}
mkdir -p "$OUT"
cd "$ROOT" || exit 2
if [ -n "${TARAL:-}" ]; then
  ENGINE=$TARAL
else
  ENGINE=$OUT/taral
  echo "== building engine (g++ -O3 -march=native -std=c++17) from src/ at $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  g++ -O3 -march=native -std=c++17 -o "$ENGINE" src/*.cpp || { echo "build failed"; exit 2; }
fi
python3 -c "import numpy, highspy" 2>/dev/null || { echo "python3 needs numpy and highspy"; exit 2; }
LIMIT=0
[ -n "${QUICK:-}" ] && LIMIT=10
rc=0
T0=$(date +%s)
declare -A HRC
for spec in "simplex all+big" "dual lp+biglp" "ipm lp+biglp"; do
  set -- $spec
  m=$1; kind=$2
  d=$OUT/$m
  mkdir -p "$d"
  echo "== harness: --method $m --kind $kind (jobs=$JOBS, cap ${TL}s/case)"
  s=$(date +%s)
  python3 tools/adv/harness.py --engine "$ENGINE" --out "$d" --jobs "$JOBS" --time-limit "$TL" --kind "$kind" --method "$m" --limit "$LIMIT" 2> "$d/harness.log"
  HRC[$m]=$?
  echo "   harness exit status ${HRC[$m]}, $(( $(date +%s) - s )) s"
  [ "${HRC[$m]}" -ne 0 ] && { echo "   harness failed for $m (see $d/harness.log)"; exit 2; }
done

echo
echo "== simplex path vs known failures (tests/known_failures.tsv)"
python3 tools/adv/report.py --results "$OUT/simplex/results.jsonl" --known tests/known_failures.tsv --md "$OUT/simplex/table.md"
r=$?; echo "   report exit status $r"; [ $r -ne 0 ] && rc=1

echo
echo "== per-category before/after (before = simplex; after = dual, ipm)"
python3 tools/adv/compare.py --before simplex="$OUT/simplex/results.jsonl" --after dual="$OUT/dual/results.jsonl" ipm="$OUT/ipm/results.jsonl" --md "$OUT/compare.md"
r=$?; echo "   compare exit status $r (1 = path-specific NEW failures)"; [ $r -ne 0 ] && rc=1

echo
echo "== reproducers (tests/repro, default path)"
python3 tools/adv/repro.py --engine "$ENGINE" --dir tests/repro
r=$?; echo "   repro exit status $r"; [ $r -ne 0 ] && rc=1

echo
echo "== method-specific reproducers (tests/repro_methods, each replayed with its own --method)"
python3 tools/adv/repro.py --engine "$ENGINE" --dir tests/repro_methods
r=$?; echo "   repro_methods exit status $r"; [ $r -ne 0 ] && rc=1

echo
echo "suite wall time $(( $(date +%s) - T0 )) s; exit status $rc (0 = no NEW failures)"
exit $rc
