#!/usr/bin/env bash
# End-to-end adversarial suite: build the engine, run the differential harness over every generated case, the MPS
# parser strictness table and the reproducer replay, then print a pass/fail table by category.
#   tests/run_adversarial.sh                 full suite (about 15 minutes on 4 CPU cores, 60 s cap per case)
#   QUICK=1 tests/run_adversarial.sh         10 cases per category (about a minute)
#   TARAL=/path/to/engine JOBS=8 OUT=dir tests/run_adversarial.sh
# Needs: g++ (C++17), python3 with numpy and highspy (HiGHS is the reference only; the engine never links it).
# Exit status: 0 when every failure is a KNOWN one (tests/known_failures.tsv, each with a reproducer), 1 on any NEW
# failure or reproducer regression, 2 when the suite could not run.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$ROOT/out/adversarial}
JOBS=${JOBS:-4}
TL=${TL:-60}
mkdir -p "$OUT"
cd "$ROOT" || exit 2

if [ -n "${TARAL:-}" ]; then
  ENGINE=$TARAL
else
  ENGINE=$OUT/taral
  echo "== building engine (g++ -O3 -march=native -std=c++17) from src/ at $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o "$ENGINE" src/*.cpp || { echo "build failed"; exit 2; }
fi
python3 -c "import numpy, highspy" 2>/dev/null || { echo "python3 needs numpy and highspy (pip install numpy highspy)"; exit 2; }
LIMIT=0
[ -n "${QUICK:-}" ] && LIMIT=10
rc=0
T0=$(date +%s)

echo "== differential harness (engine vs HiGHS $(python3 -c 'import highspy;print(highspy.Highs().version())'), jobs=$JOBS, cap ${TL}s/case)"
python3 tools/adv/harness.py --engine "$ENGINE" --out "$OUT" --jobs "$JOBS" --time-limit "$TL" --limit "$LIMIT" || { echo "harness failed"; exit 2; }
echo
echo "== results by category"
python3 tools/adv/report.py --results "$OUT/results.jsonl" --known tests/known_failures.tsv --md "$OUT/table.md" || rc=1

echo
echo "== reproducers (tests/repro)"
python3 tools/adv/repro.py --engine "$ENGINE" --dir tests/repro || rc=1

echo
echo "suite wall time $(( $(date +%s) - T0 )) s; exit status $rc (0 = only known failures)"
exit $rc
