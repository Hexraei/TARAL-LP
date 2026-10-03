#!/usr/bin/env bash
# Replays the three ill-scaled QPs against expected.tsv (status, objective to 1e-8 relative, iteration count).
# Usage: tests/ipm_numerics/check.sh ENGINE
set -u
cd "$(dirname "$0")" || exit 2
ENGINE=${1:?engine binary}
rc=0
while IFS=$'\t' read -r f st obj it; do
  case "$f" in \#*|"") continue;; esac
  read -r _ gst _ gobj _ git _ < <("$ENGINE" "$f" --method ipm)
  if python3 -c "import sys; a,b=float('$gobj'),float('$obj'); sys.exit(0 if abs(a-b)<=1e-8*max(1,abs(b)) else 1)" && [ "$gst" = "$st" ] && [ "$git" = "$it" ]; then
    echo "ok    $f $gst $gobj $git"
  else
    echo "FAIL  $f got $gst $gobj $git want $st $obj $it"; rc=1
  fi
done < expected.tsv
exit $rc
