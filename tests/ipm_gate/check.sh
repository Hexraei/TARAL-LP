#!/usr/bin/env bash
# Replays tests/ipm_gate/*.mps with --method ipm and compares the status with expected.tsv. Usage: tests/ipm_gate/check.sh ENGINE
# Model e is infeasible and has an improving direction: it must stay dual_infeasible (never unbounded).
set -u
cd "$(dirname "$0")" || exit 2
ENGINE=${1:?engine binary}
rc=0
while IFS=$'\t' read -r f want; do
  case "$f" in \#*|"") continue;; esac
  got=$("$ENGINE" "$f" --method ipm | awk 'NR==1{print $2}')
  if [ "$got" = "$want" ]; then echo "ok    $f $got"; else echo "FAIL  $f got $got want $want"; rc=1; fi
done < expected.tsv
exit $rc
