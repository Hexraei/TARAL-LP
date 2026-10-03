#!/usr/bin/env python3
"""For every MPS given (or every *.mps under a directory) zero the objective and solve the constraints alone with HiGHS
(presolve off). Prints the model status per file; 'Optimal' means the constraint set is feasible.
  python3 tools/diag_qp_ipm/feasible_check.py DIR_OR_FILES..."""
import collections, glob, os, sys
import highspy
import numpy as np

files = []
for a in sys.argv[1:]:
    files += sorted(glob.glob(os.path.join(a, "**", "*.mps"), recursive=True)) if os.path.isdir(a) else [a]
tab = collections.Counter()
for p in files:
    h = highspy.Highs(); h.setOptionValue("output_flag", False); h.setOptionValue("presolve", "off")
    h.readModel(p)
    n = h.getNumCol()
    h.changeColsCost(n, np.arange(n, dtype=np.int32), np.zeros(n))
    h.run()
    st = h.modelStatusToString(h.getModelStatus())
    tab[st] += 1
    print(os.path.basename(p), st)
print(dict(tab))
