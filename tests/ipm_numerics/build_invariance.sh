#!/usr/bin/env bash
# The interior-point path must not depend on compiler or flags. Builds src/ with several configurations and requires
# identical status / objective / iteration count for the three models (and, with FULL=1, for every model in tests/ipm_gate).
# Usage: tests/ipm_numerics/build_invariance.sh            (needs g++; clang++ is used when present)
set -u
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
T=$(mktemp -d)
rc=0; ref=""
build() {  # name compiler flags...
  local name=$1 cxx=$2; shift 2
  command -v "$cxx" >/dev/null || { echo "skip  $name ($cxx not found)"; return; }
  "$cxx" "$@" -std=c++17 -o "$T/$name" "$ROOT"/src/*.cpp 2>/dev/null || { echo "FAIL  $name did not build"; rc=1; return; }
  out=$(for f in "$ROOT"/tests/ipm_numerics/bounds_huge-*.mps; do "$T/$name" "$f" --method ipm | awk '{print $2, $4, $6}'; done)
  if [ -z "$ref" ]; then ref=$out; fi
  if [ "$out" = "$ref" ]; then echo "ok    $name"; else echo "FAIL  $name differs:"; echo "$out"; rc=1; fi
}
build gcc_O3_native g++ -O3 -march=native
build gcc_O2_native g++ -O2 -march=native
build gcc_O3 g++ -O3
build clang_O3_native clang++ -O3 -march=native
echo "reference output:"; echo "$ref"
rm -rf "$T"
exit $rc
