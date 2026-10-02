#!/usr/bin/env bash
# QP adversarial suite: build the engine, run every generated QP case (tools/advqp) against HiGHS QP, constructed
# optima and independent KKT certificates, then print the per-category table. Tests and tools only, no src/ changes.
#   tests/run_adversarial_qp.sh                 full suite (about 15 minutes on 4 cores)
#   QUICK=1 tests/run_adversarial_qp.sh         8 cases per category
#   TARAL=/path/to/engine JOBS=8 OUT=dir tests/run_adversarial_qp.sh
# Needs: g++ (C++17), python3 with numpy and highspy. Exit status: 0 when every failure is a KNOWN one
# (tests/qp_known_failures.tsv), 1 on a NEW failure or reproducer regression, 2 when the suite could not run.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$ROOT/out/adversarial_qp}
JOBS=${JOBS:-4}
TL=${TL:-30}
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
[ -n "${QUICK:-}" ] && LIMIT=8
rc=0
echo "== generator / writer / certificate self-check (no engine)"
python3 tools/advqp/qpharness.py --selftest --limit 6 || echo "   (self-check differences above are reader/reference notes, see docs/adversarial_qp.md)"
echo "== differential harness (engine vs HiGHS $(python3 -c 'import highspy;print(highspy.Highs().version())'), jobs=$JOBS, cap ${TL}s/case)"
python3 tools/advqp/qpharness.py --engine "$ENGINE" --out "$OUT" --jobs "$JOBS" --time-limit "$TL" --limit "$LIMIT" || { echo "harness failed"; exit 2; }
echo
echo "== results by category"
python3 tools/advqp/qpreport.py --results "$OUT/results.jsonl" --known tests/qp_known_failures.tsv --md "$OUT/table.md" || rc=1
if [ -d tests/qp_repro ]; then
  echo
  echo "== reproducers (tests/qp_repro)"
  python3 tools/advqp/qprepro.py --engine "$ENGINE" --dir tests/qp_repro || rc=1
fi
exit $rc
