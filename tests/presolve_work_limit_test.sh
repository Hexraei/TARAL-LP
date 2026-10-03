#!/bin/sh
set -eu
bin=${1:-build/taral}
fixture=tests/fixtures/presolve_fixedbig.mps
for flags in '--presolve --work-limit 3' '--work-limit 3 --presolve' '--presolve-log /no-such-dir/combined.json --work-limit 3'; do
    set +e
    output=$($bin "$fixture" $flags 2>&1)
    code=$?
    set -e
    [ "$code" -eq 2 ] || { echo "FAIL combined flags exit $code"; exit 1; }
    echo "$output" | grep -q 'cannot be combined'
done
echo 'PASS presolve/work-limit combinations rejected before execution'
