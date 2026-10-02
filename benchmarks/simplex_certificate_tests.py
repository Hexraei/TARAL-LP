#!/usr/bin/env python3
"""Deterministic certificate regression suite; optional SciPy random LP reference.
Usage: python benchmarks/simplex_certificate_tests.py --engine ./taral
"""
import argparse
import copy
import json
import random
import subprocess
import tempfile
from pathlib import Path
from mps_semantics_check import CASES
from simplex_certificate_check import verify

# random_13 from the reviewer: row activities land a few ulps either side of range endpoints
ENDPOINT_TIE_MPS = 'NAME RAND\nOBJSENSE MAX\nROWS\n N obj\n L r0\n L r1\n L r2\n L r3\n L r4\n L r5\n L r6\nCOLUMNS\n x0 obj 3.0\n x0 r0 2.0\n x0 r1 4.0\n x0 r2 -1.0\n x0 r3 -3.0\n x0 r4 -2.0\n x0 r5 -4.0\n x0 r6 -5.0\n x1 obj -1.0\n x1 r0 -4.0\n x1 r1 3.0\n x1 r2 -2.0\n x1 r3 -3.0\n x1 r4 -1.0\n x1 r5 -3.0\n x1 r6 -5.0\n x2 obj 4.0\n x2 r1 1.0\n x2 r2 -5.0\n x2 r3 -2.0\n x2 r4 2.0\n x2 r5 -2.0\n x2 r6 -5.0\n x3 obj 4.0\n x3 r1 -1.0\n x3 r2 1.0\n x3 r5 -3.0\n x3 r6 -1.0\n x4 obj 4.0\n x4 r0 4.0\n x4 r1 -3.0\n x4 r3 -4.0\n x4 r4 5.0\n x4 r5 1.0\n x4 r6 -4.0\n x5 obj -4.0\n x5 r0 2.0\n x5 r1 -4.0\n x5 r2 -5.0\n x5 r3 5.0\n x5 r4 -4.0\n x5 r5 2.0\n x5 r6 -2.0\nRHS\n rhs obj 9\n rhs r0 12.0\n rhs r1 -1.0\n rhs r2 -2.5\n rhs r3 6.0\n rhs r4 19.0\n rhs r5 9.0\n rhs r6 7.5\nRANGES\n rng r0 9.0\n rng r1 15.0\n rng r2 5.0\n rng r3 10.0\n rng r4 5.0\n rng r5 8.0\n rng r6 13.0\nBOUNDS\n LO b x0 -4.0\n UP b x0 2.0\n LO b x1 -3.0\n UP b x1 -2.0\n LO b x2 0.0\n UP b x2 6.0\n LO b x3 -2.0\n UP b x3 2.0\n LO b x4 0.0\n UP b x4 1.0\n LO b x5 -4.0\n UP b x5 3.0\nENDATA\n'

EXTRA = {
    'endpoint_roundoff_tie': ENDPOINT_TIE_MPS,
    'no_rows_boxed': 'NAME BOX\nROWS\n N obj\nCOLUMNS\n x obj -2\n y obj 3\nRHS\n rhs obj -4\nBOUNDS\n UP b x 5\n LO b y -2\nENDATA\n',
    'no_columns': 'NAME EMPTY\nROWS\n N obj\n L r\nRHS\n rhs r 1 obj -5\nENDATA\n',
    'zero_cost_free': 'NAME FREE\nROWS\n N obj\nCOLUMNS\n x obj 0\nBOUNDS\n FR b x\nENDATA\n',
    'fixed_column_dual': 'NAME FIX\nROWS\n N obj\n E r\nCOLUMNS\n x obj -3 r 1\n y obj 2 r 1\nRHS\n rhs r 4\nBOUNDS\n FX b x 1\nENDATA\n',
    'degenerate': 'NAME DEG\nROWS\n N obj\n G r1\n G r2\nCOLUMNS\n x obj 1 r1 1\n x r2 2\nRHS\n rhs r1 1 r2 2\nENDATA\n',
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--engine', default='./taral')
    a = ap.parse_args()
    total = rejected = 0
    with tempfile.TemporaryDirectory() as directory:
        d = Path(directory)
        def run(name, text, reference=None):
            nonlocal total, rejected
            mps, js = d / (name+'.mps'), d / (name+'.json')
            mps.write_text(text)
            for method in ('simplex', 'dual'):
                subprocess.run([a.engine, str(mps), '--method', method, '--json', str(js)], check=True, capture_output=True, timeout=15)
                cert = json.loads(js.read_text())
                result = verify(mps, cert)
                assert cert['certificate_quality']=='pass'
                assert result['pass_certificate'], (name, method, cert, result)
                if reference is not None:
                    assert abs(cert['objective']-reference) <= 1e-8*(1+abs(reference)), (name, method, reference, cert)
                total += 1
                # An independent checker must notice edits even if residual fields stay zero.
                for key in ('objective', 'dual_objective', 'primal_res', 'dual_res', 'gap', 'complementarity'):
                    bad = copy.deepcopy(cert); bad[key] += 1
                    assert not verify(mps, bad)['pass_certificate'], key
                    rejected += 1
                for key in ('x', 'row_dual', 'reduced_cost', 'row_activity', 'row_violation_abs', 'row_term_magnitude', 'row_violation_magnitude_scaled'):
                    if cert[key]:
                        bad = copy.deepcopy(cert); bad[key][0] += 1
                        valid_alternate = name == 'zero_cost_free' and key == 'x'
                        if not valid_alternate:
                            assert not verify(mps, bad)['pass_certificate'], (name, key)
                            rejected += 1
                        bad[key][0] = None
                        try:
                            passed = verify(mps, bad)['pass_certificate']
                        except ValueError:
                            passed = False
                        assert not passed
                        rejected += 1
            print('PASS', name)
        for name, (text, status, objective, _) in CASES.items():
            if status == 'optimal' and not any(v in text for v in ("'MARKER'", 'QUADOBJ', 'QMATRIX', 'QSECTION', '\n BV ', '\n LI ', '\n UI ')):
                run(name, text, objective)
        for name, text in EXTRA.items(): run(name, text)
        try:
            from scipy.optimize import linprog
            import numpy as np
        except ImportError:
            print('SKIP random SciPy reference (SciPy/NumPy not installed)')
        else:
            rng = random.Random(57463)
            for test in range(80):
                n, m = rng.randrange(2, 9), rng.randrange(1, 9)
                A = np.array([[rng.randint(-5, 5) for _ in range(n)] for _ in range(m)], dtype=float)
                c = np.array([rng.randint(-4, 4) for _ in range(n)], dtype=float)
                lo = np.array([rng.randint(-4, 0) for _ in range(n)], dtype=float)
                up = lo + np.array([rng.randint(1, 7) for _ in range(n)])
                witness = (lo+up)/2
                center = A @ witness
                lower = center - np.array([rng.randint(0, 8) for _ in range(m)])
                upper = center + np.array([rng.randint(0, 8) for _ in range(m)])
                sense = -1 if test % 2 else 1
                ref = linprog(sense*c, A_ub=np.vstack([A, -A]), b_ub=np.r_[upper, -lower], bounds=list(zip(lo,up)), method='highs')
                assert ref.success, ref.message
                offset = rng.randint(-9, 9)
                lines = ['NAME RAND', 'OBJSENSE '+('MAX' if sense < 0 else 'MIN'), 'ROWS', ' N obj']
                lines += [f' L r{i}' for i in range(m)]
                lines += ['COLUMNS']
                for j in range(n):
                    lines += [f' x{j} obj {c[j]}'] + [f' x{j} r{i} {A[i,j]}' for i in range(m) if A[i,j]]
                lines += ['RHS', f' rhs obj {-offset}']+[f' rhs r{i} {upper[i]}' for i in range(m)]
                lines += ['RANGES']+[f' rng r{i} {upper[i]-lower[i]}' for i in range(m)]
                lines += ['BOUNDS']
                for j in range(n): lines += [f' LO b x{j} {lo[j]}', f' UP b x{j} {up[j]}']
                lines += ['ENDATA']
                run('random_'+str(test), '\n'.join(lines)+'\n', sense*ref.fun+offset)
    print(f'PASS {total} simplex/dual certificates; {rejected} corrupted/nonfinite certificates rejected')


if __name__ == '__main__':
    main()
