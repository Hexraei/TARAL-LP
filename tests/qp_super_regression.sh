#!/usr/bin/env bash
# Regression for the opt-in supernodal QP kernel (TARAL_QP_SUPER=1, added October 16, 2026):
#   1. knob parsing must not perturb the default path: TARAL_QP_SUPER=1 at the default 2e8-flop
#      cutoff (supernodal stays off on small models) must match the unset run on every field but wall time.
#   2. forced supernodal (TARAL_QP_SUPER_MIN_FLOPS=0) must return the same status and objective as the
#      scalar path on a tiny convex QP with a known analytic optimum (min x1^2+x2^2 s.t. x1+x2>=2, xi<=10; obj 2).
# The scalar fallback on a non-finite supernodal factor is covered by the audited run record in
# results/qplib_qp_super_20261016/ (QPLIB_10038), not by this quick test.
#   tests/qp_super_regression.sh            build and run (under a minute)
#   TARAL=/path/to/engine tests/qp_super_regression.sh
# Exit 0 pass, 1 mismatch, 2 could not run.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$ROOT/out/qp_super_regression}
mkdir -p "$OUT"
cd "$ROOT" || exit 2
if [ -n "${TARAL:-}" ]; then
  ENGINE=$TARAL
else
  ENGINE=$OUT/taral
  echo "== building engine (g++ -O2 -std=c++17 -pthread) from src/ at $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  g++ -O2 -std=c++17 -pthread -o "$ENGINE" src/*.cpp || { echo "build failed"; exit 2; }
fi
cat > "$OUT/tiny_qp.mps" <<'MPS'
NAME          TINYQP
ROWS
 N  obj
 G  c1
COLUMNS
    x1        c1        1.0
    x2        c1        1.0
RHS
    rhs       c1        2.0
BOUNDS
 UP bnd       x1        10.0
 UP bnd       x2        10.0
QUADOBJ
    x1        x1        2.0
    x2        x2        2.0
ENDATA
MPS
"$ENGINE" "$OUT/tiny_qp.mps" --method ipm --json "$OUT/unset.json" >/dev/null 2>&1 || { echo "unset run failed"; exit 2; }
TARAL_QP_SUPER=1 "$ENGINE" "$OUT/tiny_qp.mps" --method ipm --json "$OUT/knob.json" >/dev/null 2>&1 || { echo "knob run failed"; exit 2; }
TARAL_QP_SUPER=1 TARAL_QP_SUPER_MIN_FLOPS=0 "$ENGINE" "$OUT/tiny_qp.mps" --method ipm --json "$OUT/forced.json" >/dev/null 2>&1 || { echo "forced run failed"; exit 2; }
python3 - "$OUT" <<'PY'
import json, sys
out = sys.argv[1]
u = json.load(open(f"{out}/unset.json"))
k = json.load(open(f"{out}/knob.json"))
f = json.load(open(f"{out}/forced.json"))
fields = ["status", "objective", "dual_objective", "primal_res", "dual_res", "gap",
          "max_row_viol", "max_bound_viol", "iterations", "message"]
rc = 0
for fld in fields:
    if u.get(fld) != k.get(fld):
        print(f"MISMATCH knob-unset identity field {fld}: {u.get(fld)} vs {k.get(fld)}"); rc = 1
if abs(u["objective"] - 2.0) > 1e-6:
    print(f"MISMATCH scalar objective vs analytic optimum 2.0: {u['objective']}"); rc = 1
if f["status"] != u["status"]:
    print(f"MISMATCH forced-super status: {f['status']} vs {u['status']}"); rc = 1
if f["objective"] is None or abs(f["objective"] - u["objective"]) > 1e-9 * (1 + abs(u["objective"])):
    print(f"MISMATCH forced-super objective: {f['objective']} vs {u['objective']}"); rc = 1
print("qp_super_regression:", "PASS" if rc == 0 else "FAIL")
sys.exit(rc)
PY
