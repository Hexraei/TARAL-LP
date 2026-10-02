#!/usr/bin/env bash
# Independent second-host run. Prints everything to stdout (Render logs are the only output channel); ends with SECOND_HOST_DONE then idles so logs can be saved. Delete the service after saving logs.
set -u; cd /app
echo "SECOND_HOST_BEGIN $(date -u +%FT%TZ)"; echo "== host"; uname -a; lscpu | grep -E 'Model name|^CPU\(s\)|^Thread|MHz' ; free -g | head -2; g++ --version | head -1; python3 -c "import numpy,scipy;print('numpy',numpy.__version__,'scipy',scipy.__version__)"
echo "== package sha256"; sha256sum src/* | sort -k2
echo "BUILD g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp"
g++ -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic -o taral src/*.cpp 2>&1 | tail -20; echo "BUILD_RC=${PIPESTATUS[0]}"
# 1. 93-case netlib gate, pinned corpus, 60 s cap (public sanitized gate)
bash gate/fetch_corpus.sh /app/corpus
for MODE in primal dual; do
  EA=""; [ $MODE = dual ] && EA="--method dual"
  echo "== GATE $MODE"; python3 gate/gate_pub.py --cmd /app/taral --corpus /app/corpus --cases "$(cat gate/cases93.txt)" --time-limit 60 --label host2_$MODE --out /app/out_$MODE --engine-args "$EA" 2>&1 | tail -40
done
# 2. adversarial suite (from the test tree shipped in tests/)
if [ -f tests/run_adversarial.sh ]; then echo "== ADVERSARIAL"; (bash tests/run_adversarial.sh 2>&1 | tail -80; echo "ADVERSARIAL_EXIT=${PIPESTATUS[0]}"; cat out/adversarial/table.md out/adversarial/parser_table.md 2>/dev/null | head -150); fi
# 3. Mittelmann big/memory cases, solver only (no HiGHS reference here: 600 s HiGHS on 4 cores is the Kaggle protocol), 300 s cap
if [ -f mit/emps.c ]; then gcc -O2 -o /app/emps mit/emps.c; echo "== MIT"; python3 mit/mit_run.py /app H2A 300 0 2>&1 | tail -60; fi
echo "SECOND_HOST_DONE $(date -u +%FT%TZ)"; sleep infinity
