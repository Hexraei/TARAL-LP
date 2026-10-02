#!/usr/bin/env bash
# One command: build the C++ engine, fetch the pinned Netlib corpus, run the gate, print figures.
# Usage: bash reproduce.sh [cases comma list]   (default: all 93 + truss, 60 s cap per case)
# Env: CORPUS=<dir of .mps> to skip the clone; TL=60 time cap; OUT=out
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT/reproduction"
CANONICAL="$ROOT/cpp-engine"
TL=${TL:-60}; OUT=${OUT:-out}
mkdir -p "$OUT"
REPO=https://github.com/ozy4dm/lp-data-netlib.git
PIN=56257eea85b433ce6aa67d26156b36385318fd6f
echo "== environment"; python3 --version; g++ --version | head -1
python3 - <<'VERSIONS'
from importlib.metadata import version, PackageNotFoundError
for name in ('numpy', 'scipy', 'highspy'):
    try: print(name, version(name))
    except PackageNotFoundError: print(name, 'not installed')
VERSIONS
echo "== cold build (reported separately from solver time)"
s=$(date +%s.%N); g++ -O2 -std=c++17 -o "$OUT/taral" "$CANONICAL/engine/taral.cpp"; e=$(date +%s.%N)
echo "compile_seconds=$(python3 -c "print(round($e-$s,1))")"
echo "engine_sha256=$(sha256sum "$CANONICAL/engine/taral.cpp" | cut -d' ' -f1)"
echo "backend=CPU (single-thread C++ simplex; no GPU backend in this engine)"
if [ -z "${CORPUS:-}" ]; then
  echo "== fetching pinned corpus ($PIN)"
  rm -rf "$OUT/netlib_corpus" && git clone -q "$REPO" "$OUT/netlib_corpus" && git -C "$OUT/netlib_corpus" checkout -q "$PIN"
  CORPUS="$OUT/netlib_corpus/mps_files"
fi
CORPUS=$(cd "$CORPUS" && pwd)
MANIFEST="$ROOT/reproduction/MANIFEST_corpus_sha256.txt"
echo "== corpus file hashes vs manifest"
(cd "$CORPUS" && sha256sum -c --quiet "$MANIFEST" 2>&1 | tail -3)
echo "manifest check done"
echo "== gate (protocol ${TL}s)"
CASES=()
[ $# -ge 1 ] && CASES=(--cases "$1") || CASES=(--with-extra)
python3 "$CANONICAL/tools/gate_cpp.py" --cmd "$OUT/taral" --corpus "$CORPUS" --time-limit "$TL" --protocol "${TL}s" --out "$OUT" "${CASES[@]}"
