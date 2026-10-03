#!/usr/bin/env python3
"""Classify the QPLIB instance table. Source: https://qplib.zib.de/instances.html (columns in the table header:
Cvx = continuous relaxation convex, O = objective type, V = variable type, C = constraint type, Quad.Cons. = number of quadratic constraints).
Rule, applied in this order (exclusive classes):
  in_class            : Cvx checked AND O in {C,D} (convex / diagonal-convex quadratic objective) AND V == C (all continuous)
                        AND C in {L,B,N} (linear, box or no constraints) AND Quad.Cons. == 0
  quadratic_constraints: Quad.Cons. > 0
  integer_binary      : V != C (B binary, M mixed, I integer, G general)
  nonconvex_continuous: remaining (not Cvx)
Usage: python3 qplib_classify.py [instances.html] > qplib_classification.csv"""
import re, sys, csv, urllib.request, collections
s = open(sys.argv[1]).read() if len(sys.argv) > 1 else urllib.request.urlopen("https://qplib.zib.de/instances.html").read().decode()
rows = re.findall(r"<TR[^>]*>(.*?)</TR>", s[s.find("<TBODY>"):], re.S)
w = csv.writer(sys.stdout); w.writerow(["qplib_id", "convex", "objective_type", "variable_type", "constraint_type", "n_vars", "n_cons", "n_quad_cons", "class"])
cnt = collections.Counter()
for r in rows:
    td = [re.sub(r"<[^>]+>", "", x).replace("&nbsp;", "").strip() for x in re.findall(r"<TD[^>]*>(.*?)</TD>", r, re.S)]
    m = re.search(r"QPLIB_(\d+)\.html", r)
    if not m or len(td) < 13:
        continue
    cvx, O, V, C, nv, nc, nq = "10004" in td[1], td[2], td[5], td[9], td[6], td[10], td[11]
    if cvx and O in ("C", "D") and V == "C" and C in ("L", "B", "N") and nq in ("0", ""):
        k = "in_class"
    elif nq not in ("0", ""):
        k = "quadratic_constraints"
    elif V != "C":
        k = "integer_binary"
    else:
        k = "nonconvex_continuous"
    cnt[k] += 1
    w.writerow([m.group(1), cvx, O, V, C, nv, nc, nq or 0, k])
print("# total", sum(cnt.values()), dict(cnt), file=sys.stderr)
