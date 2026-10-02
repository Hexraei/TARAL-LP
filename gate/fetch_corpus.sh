#!/usr/bin/env bash
# fetch pinned netlib corpus and verify sha256 manifest
set -u; D=${1:-/app/corpus}; mkdir -p $D; SHA=56257eea85b433ce6aa67d26156b36385318fd6f; bad=0
while read -r h n; do n=${n#\*}; curl -sfL --retry 3 -o $D/$n https://raw.githubusercontent.com/ozy4dm/lp-data-netlib/$SHA/mps_files/$n || echo "FETCH_FAIL $n"
  [ "$(sha256sum $D/$n | cut -d' ' -f1)" = "$h" ] || { echo "HASH_MISMATCH $n"; bad=$((bad+1)); }; done < "$(dirname "$0")/manifest-sha256.txt"
echo "corpus hash mismatches: $bad"
