#!/usr/bin/env bash
# Run in a Codespace from the branch root:  nohup bash codespace_run.sh > run.log 2>&1 &
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo -n"
$SUDO apt-get update -qq
$SUDO apt-get install -y -qq --no-install-recommends g++ gcc make python3 python3-pip python3-venv curl ca-certificates bzip2 gzip git util-linux
python3 -m venv "$HOME/venv"
"$HOME/venv/bin/pip" install -q --no-cache-dir numpy scipy==1.15.3 highspy==1.15.1
export PATH="$HOME/venv/bin:$PATH"
python3 -c "import scipy,highspy;assert scipy.__version__=='1.15.3',scipy.__version__"
$SUDO rm -rf /app; $SUDO ln -s "$PWD" /app
# fail loud: make build and corpus fetch fatal in a working copy of run_all.sh
sed -e 's/^set -u; cd \/app/set -u; set -o pipefail; cd \/app/' \
    -e 's#^g++ \(.*\) 2>&1 | tail -20; echo "BUILD_RC=\${PIPESTATUS\[0\]}"#g++ \1 2>\&1 | tail -20; BRC=${PIPESTATUS[0]}; echo "BUILD_RC=$BRC"; [ "$BRC" = 0 ] || { echo RUN_FAILED_BUILD; exit 1; }#' \
    -e 's#^bash gate/fetch_corpus.sh /app/corpus#bash gate/fetch_corpus.sh /app/corpus || { echo RUN_FAILED_CORPUS; exit 1; }#' \
    -e 's/sleep infinity//' run_all.sh > /tmp/run_all_cs.sh
grep -q RUN_FAILED_BUILD /tmp/run_all_cs.sh && grep -q RUN_FAILED_CORPUS /tmp/run_all_cs.sh || { echo "PATCH_FAILED"; exit 1; }
bash /tmp/run_all_cs.sh 2>&1 | tee full_run.log
grep -q '^SECOND_HOST_DONE' full_run.log || { echo RUN_FAILED_INCOMPLETE; exit 1; }
mkdir -p results; cp full_run.log results/; cp -r out_primal out_dual out results/ 2>/dev/null || true
git config user.name >/dev/null && git config user.email >/dev/null || { echo RUN_FAILED_NO_GIT_IDENTITY; exit 1; }
git add results
git commit -q -m "second-host results"
git push -q origin HEAD
echo COMMITTED_RESULTS
