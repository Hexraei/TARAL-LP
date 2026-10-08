#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
output=${1:-/tmp/r2-miqp}
c++ -O2 -std=c++17 -I"$root/src" -I"$root/experiments/r2_miqp" \
 "$root/experiments/r2_miqp/main.cpp" "$root/experiments/r2_miqp/miqp.cpp" \
 "$root/src/ipm.cpp" "$root/src/mps.cpp" -o "$output"
