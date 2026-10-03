#!/bin/sh
# CLI validation tests. Usage: tests/cli_test.sh [path/to/taral]   (default: ./taral)
BIN=${1:-./taral}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
fail=0
cat > "$T/lp.mps" <<'M'
NAME          T
ROWS
 N  OBJ
 L  C1
COLUMNS
    X         OBJ       -1.0   C1        1.0
RHS
    RHS       C1        4.0
ENDATA
M
# expect NAME CODE [args...]
expect() { name=$1; code=$2; shift 2
  "$BIN" "$@" >"$T/out" 2>"$T/err"; rc=$?
  if [ "$rc" = "$code" ]; then echo "PASS $name (exit $rc)"; else echo "FAIL $name: exit $rc, want $code"; fail=1; fi; }
expect help 0 --help
grep -q '^usage: taral' "$T/out" || { echo "FAIL help text"; fail=1; }
expect help-short 0 -h
expect no-args 2
expect valid-default 0 "$T/lp.mps"
grep -q 'status optimal' "$T/out" || { echo "FAIL valid-default status"; fail=1; }
expect valid-opts 0 "$T/lp.mps" --time-limit 5 --node-limit 10 --method dual --sol "$T/s" --json "$T/j"
for m in simplex dual ipm; do expect "method-$m" 0 "$T/lp.mps" --method $m; grep -q 'status optimal' "$T/out" || { echo "FAIL method-$m status"; fail=1; }; done
expect bad-method 2 "$T/lp.mps" --method bogus
expect missing-model-file 3 /nonexistent.mps
grep -q parse_error "$T/out" || { echo "FAIL missing file keeps parse_error"; fail=1; }
for v in 0 -1 abc 1x nan inf -inf 1e999 ""; do expect "time-limit[$v]" 2 "$T/lp.mps" --time-limit "$v"; done
for v in -1 abc 1.5 "" 99999999999999999999999; do expect "node-limit[$v]" 2 "$T/lp.mps" --node-limit "$v"; done
expect time-limit-missing 2 "$T/lp.mps" --time-limit
expect sol-missing 2 "$T/lp.mps" --sol
expect json-missing 2 "$T/lp.mps" --json
expect unknown-option 2 "$T/lp.mps" --bogus
expect extra-arg 2 "$T/lp.mps" other.mps
expect double-dash 0 -- "$T/lp.mps"
expect help-after-model 0 "$T/lp.mps" --help
exit $fail
