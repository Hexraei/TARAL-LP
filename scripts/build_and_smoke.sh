#!/usr/bin/env bash
# CPU-only build and small, deterministic checks. No Netlib download or CUDA needed.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/build}"
PYTHON="${PYTHON:-python3}"
cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --parallel "${BUILD_JOBS:-2}"
"$PYTHON" -m venv "$BUILD_DIR/smoke-venv"
"$BUILD_DIR/smoke-venv/bin/python" -m pip install -r "$ROOT/benchmarks/requirements-smoke.txt"
"$BUILD_DIR/smoke-venv/bin/python" "$ROOT/benchmarks/smoke.py" \
    --engine "$BUILD_DIR/taral" --out "$BUILD_DIR/smoke" "$@"
