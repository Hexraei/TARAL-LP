#!/usr/bin/env python3
"""Five equality-space labels must certify; negative reduced curvature must not pass.
Usage: python tests/ipm_gate/check_equality_labels.py ENGINE
Needs numpy and highspy, like tools/advqp/qpharness.py.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/advqp'))
import qpcase
import qpequality
import qpharness
with tempfile.TemporaryDirectory() as temp:
    for i in (12, 27, 34, 60, 76):
        case = qpcase.gen('indefinite', i)
        model, sol, js = (Path(temp) / n for n in ('case.mps', 'case.sol', 'case.json'))
        model.write_text(qpcase.write_mps(case))
        subprocess.run([sys.argv[1], str(model), '--method', 'ipm', '--sol', str(sol), '--json', str(js)], check=True, capture_output=True)
        ours = json.loads(js.read_text())
        x = qpharness.read_sol(sol, case.n)
        verdict = qpharness.judge(case, ours, x, 30)
        assert verdict['verdict'] == 'PASS', (case.id, verdict)
        case.Q = (1.0 if case.maximize else -1.0) * np.eye(case.n)
        assert qpequality.reduced_case(case) is None
        x = np.zeros(case.n)
        bad = qpharness.judge(case, {'status': 'optimal', 'objective': case.obj(x)}, x, 30)
        assert bad['verdict'] == 'FAIL_STATUS', (case.id, bad)
        print('ok', case.id, 'point/objective/KKT and negative-curvature guard')
    # No unrelated identity receives the scoped correction.
    case = qpcase.gen('indefinite', 12)
    case.id = 'indefinite-99999'
    assert qpequality.reduced_case(case) is None
