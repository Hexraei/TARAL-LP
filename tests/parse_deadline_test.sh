#!/bin/sh
# Small-input regression for parse timeout status/JSON, without a large MPS corpus.
set -eu
BIN=${1:-./taral}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cat > "$T/lp.mps" <<'M'
NAME T
ROWS
 N OBJ
 L C1
COLUMNS
 X OBJ -1 C1 1
RHS
 RHS C1 4
ENDATA
M
set +e
"$BIN" "$T/lp.mps" --time-limit 0.000000001 --json "$T/result.json" > "$T/out" 2> "$T/err"
rc=$?
set -e
test "$rc" = 4
grep -q '^status time_limit: parse exceeded the time limit' "$T/out"
python3 - "$T/result.json" <<'PY'
import json,sys
r=json.load(open(sys.argv[1]))
assert r['status']=='time_limit',r
assert r['message']=='parse exceeded the time limit',r
PY
"$BIN" "$T/lp.mps" > "$T/default"
grep -q '^status optimal' "$T/default"
echo 'PASS parse-deadline timeout exit/JSON and unlimited small-input parse'
