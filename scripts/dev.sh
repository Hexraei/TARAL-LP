#!/usr/bin/env bash
# One entry point for common developer tasks. Run from anywhere; paths resolve from the repo root.
#   scripts/dev.sh build    CMake Release build with warnings into build/ (build/taral)
#   scripts/dev.sh test     build with unit tests enabled and run ctest
#   scripts/dev.sh pinned   the exact measurement build used for the committed ledgers: ./taral
#   scripts/dev.sh smoke    CPU-only build and small deterministic checks (scripts/build_and_smoke.sh)
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/build}"
JOBS="${BUILD_JOBS:-2}"
cmd="${1:-}"; shift || true
case "$cmd" in
  build)
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DTARAL_WARNINGS=ON
    cmake --build "$BUILD_DIR" --parallel "$JOBS" ;;
  test)
    cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DTARAL_WARNINGS=ON -DTARAL_BUILD_TESTS=ON
    cmake --build "$BUILD_DIR" --parallel "$JOBS"
    ctest --test-dir "$BUILD_DIR" --output-on-failure ;;
  pinned)
    # Same command as ci/smoke.sh and reproduction/: do not change flags here.
    cd "$ROOT" && g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp ;;
  smoke)
    exec "$ROOT/scripts/build_and_smoke.sh" "$@" ;;
  *)
    sed -n '2,6p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2; exit 2 ;;
esac
