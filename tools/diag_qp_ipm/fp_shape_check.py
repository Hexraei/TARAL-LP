#!/usr/bin/env python3
"""Compare the floating-point instruction shape of two builds of src/ipm.cpp, symbol by symbol.
For every function (lambdas included) present in both object files, the sequence of FP arithmetic opcodes
(mul/add/sub/div/sqrt and every fused multiply-add form) is extracted; any symbol whose sequence differs between the
two builds can change the last bits of the interior-point iterates. Register allocation and non-FP code are ignored.
  g++ -O3 -march=native -std=c++17 -c -o a.o a/ipm.cpp ; g++ ... -o b.o b/ipm.cpp
  python3 tools/diag_qp_ipm/fp_shape_check.py a.o b.o [symbol-substring]"""
import re, subprocess, sys

FP = re.compile(r"^v?(f?n?m(add|sub)\d*[a-z]*|mul[sp]d|add[sp]d|sub[sp]d|div[sp]d|sqrt[sp]d)$")

def shapes(obj):
    out = subprocess.run(["objdump", "-d", "--no-show-raw-insn", "-C", obj], capture_output=True, text=True).stdout
    res, cur = {}, None
    for line in out.splitlines():
        m = re.match(r"^[0-9a-f]+ <(.*)>:$", line)
        if m:
            cur = m.group(1)
            res[cur] = []
            continue
        m = re.match(r"^ +[0-9a-f]+:\t(\S+)", line)
        if m and cur is not None and FP.match(m.group(1)):
            res[cur].append(m.group(1))
    return res

a, b = shapes(sys.argv[1]), shapes(sys.argv[2])
sub = sys.argv[3] if len(sys.argv) > 3 else ""
bad = n = 0
for name in sorted(set(a) & set(b)):
    if sub not in name or "cold" in name:
        continue
    n += 1
    if a[name] != b[name]:
        bad += 1
        print("DIFFERS  %3d -> %3d FP ops (%d -> %d fused): %s" % (len(a[name]), len(b[name]),
              sum("fm" in o or "fnm" in o for o in a[name]), sum("fm" in o or "fnm" in o for o in b[name]), name[:110]))
only = [k for k in set(a) ^ set(b) if sub in k and "cold" not in k and (a.get(k) or b.get(k))]
for k in sorted(only):
    print("ONLY IN %s with FP ops: %s" % ("first" if k in a else "second", k[:110]))
print("%d common symbols compared, %d differ" % (n, bad))
sys.exit(1 if bad else 0)
