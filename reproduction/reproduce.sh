#!/usr/bin/env bash
# One command: build the C++ engine, fetch the pinned Netlib corpus, run the gate, print figures.
# Usage: bash reproduce.sh [cases comma list]   (default: all 93 + truss, 60 s cap per case)
# Env: CORPUS=<dir of .mps> to skip the clone; TL=60 time cap; OUT=out
set -euo pipefail
cd "$(dirname "$0")"
TL=${TL:-60}; OUT=${OUT:-out}
REPO=https://github.com/ozy4dm/lp-data-netlib.git
PIN=56257eea85b433ce6aa67d26156b36385318fd6f
echo "== environment"; python3 --version; g++ --version | head -1
python3 -c "import numpy,scipy,highspy;print('numpy',numpy.__version__,'scipy',scipy.__version__,'highspy',getattr(highspy,'__version__','?'))" || true
echo "== cold build (reported separately from solver time)"
s=$(date +%s.%N); g++ -O2 -std=c++17 -o taral src/taral.cpp; e=$(date +%s.%N)
echo "compile_seconds=$(python3 -c "print(round($e-$s,1))")"
echo "engine_sha256=$(sha256sum src/taral.cpp | cut -d' ' -f1)"
echo "backend=CPU (single-thread C++ simplex; no GPU backend in this engine)"
if [ -z "${CORPUS:-}" ]; then
  echo "== fetching pinned corpus ($PIN)"
  rm -rf netlib_corpus && git clone -q "$REPO" netlib_corpus && git -C netlib_corpus checkout -q "$PIN"
  CORPUS=netlib_corpus
fi
echo "== corpus file hashes vs manifest"
(cd "$CORPUS" && sha256sum -c --quiet "$OLDPWD/MANIFEST_corpus_sha256.txt" 2>&1 | tail -3) && echo "manifest check done"
echo "== gate (protocol ${TL}s)"
CASES=()
[ $# -ge 1 ] && CASES=(--cases "$1") || CASES=(--with-extra)
python3 tools/gate_cpp.py --cmd ./taral --corpus "$CORPUS" --time-limit "$TL" --protocol "${TL}s" --out "$OUT" "${CASES[@]}"
