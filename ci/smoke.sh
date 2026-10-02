#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out/ci
# Exact measurement build flags; CI is correctness-only, not a timing comparison.
g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp 2>out/ci/build.log
cat out/ci/build.log
test ! -s out/ci/build.log
python ci/fetch_corpus.py
CASES=afiro,sc50a,sc50b,sc105,adlittle,blend,kb2,share2b,stocfor1,scagr7
for method in primal dual; do
  args=()
  if [ "$method" = dual ]; then args=(--engine-args='--method dual'); fi
  timeout 180 python ci/gate.py --cmd ./taral --corpus out/ci/corpus --cases "$CASES" --time-limit 10 --label "ci-$method" --out "out/ci/$method" "${args[@]}"
  # The measurement tool reports results but does not itself fail on a non-pass.
  python - "$method" <<'PY'
import csv,sys
rows=list(csv.DictReader(open('out/ci/'+sys.argv[1]+'/ledger.csv')))
assert len(rows)==10 and all(r['passed']=='True' for r in rows), [(r['case'],r['verdict']) for r in rows]
print(sys.argv[1]+': 10/10 checked optima')
PY
done
# Already-merged adversarial edge suites. Honest unresolved cases stay named in logs.
timeout 120 python benchmarks/mps_semantics_check.py --engine ./taral
timeout 120 python benchmarks/milp_tests.py --engine ./taral --per-suite 5 --seed 31 --workers 1
timeout 120 python benchmarks/qp_tests.py --engine ./taral --per-suite 3 --seed 31 --workers 1
echo 'Extended cloud adversarial suite is not merged; current coverage is parser/MILP/QP edge smoke.'
