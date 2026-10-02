"""Public gate: runs an engine that follows the CLI  ENGINE MODEL.mps --time-limit S --sol OUT.sol --json OUT.json [extra args]
on the pinned Netlib files and grades each result against HiGHS (via scipy, tooling only) solved on the ORIGINAL MPS.
Pass: engine status optimal; HiGHS status optimal; |engine obj - ref| <= max(1e-6, 1e-7|ref|); independently recomputed objective (own MPS reader, includes the
objective-row constant) agrees to the same tolerance; max relative row violation <= 1e-6; max bound violation <= 1e-6. Strict (reported only): 1e-8 on all three.
Relative row violation = violation / (1 + |violated row bound|). Denominator and tolerances are fixed. HiGHS-failed cases never count as passes."""
import sys, os, json, time, csv, subprocess, argparse
import numpy as np, scipy.sparse as sps
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpsio
from scipy.optimize import linprog
ap = argparse.ArgumentParser(); ap.add_argument('--cmd', required=True); ap.add_argument('--corpus', required=True); ap.add_argument('--cases', required=True)
ap.add_argument('--extra-case', default=''); ap.add_argument('--time-limit', type=float, default=60); ap.add_argument('--label', default='run'); ap.add_argument('--out', default='gate_out')
ap.add_argument('--fixed-highs-failed', default='pilot.we,pilot4'); ap.add_argument('--engine-args', default=''); a = ap.parse_args()
cases = [x for x in a.cases.split(',') if x]; extra = [x for x in a.extra_case.split(',') if x]; fixed = set(a.fixed_highs_failed.split(','))
os.makedirs(a.out, exist_ok=True); rows = []
def highs(A, c, rl, ru, l, u, const):
    fl = np.isfinite(rl); fh = np.isfinite(ru); eq = fl & fh & (rl == ru); ih = fh & ~eq; il = fl & ~eq
    Aub = sps.vstack([A[ih], -A[il]]).tocsr(); bub = np.concatenate([ru[ih], -rl[il]])
    r = linprog(c, A_ub=Aub if Aub.shape[0] else None, b_ub=bub if Aub.shape[0] else None, A_eq=A[eq] if eq.sum() else None, b_eq=rl[eq] if eq.sum() else None,
                bounds=list(zip(np.where(np.isfinite(l), l, None), np.where(np.isfinite(u), u, None))), method='highs', options={'time_limit': 30})
    if r.status == 2:
        r = linprog(c, A_ub=Aub if Aub.shape[0] else None, b_ub=bub if Aub.shape[0] else None, A_eq=A[eq] if eq.sum() else None, b_eq=rl[eq] if eq.sum() else None,
                    bounds=list(zip(np.where(np.isfinite(l), l, None), np.where(np.isfinite(u), u, None))), method='highs', options={'time_limit': 30, 'presolve': False})
    return ('Optimal', float(r.fun) + const) if r.status == 0 else ('Failed%d' % r.status, None)
def grade(name, in_den):
    path = os.path.join(a.corpus, name + '.mps'); sol = os.path.join(a.out, name + '.sol'); js = os.path.join(a.out, name + '.json')
    for f in (sol, js):
        if os.path.exists(f): os.remove(f)
    t = time.time(); to = False
    try: subprocess.run(a.cmd.split() + [path, '--time-limit', str(a.time_limit), '--sol', sol, '--json', js] + a.engine_args.split(), timeout=a.time_limit + 15, capture_output=True)
    except subprocess.TimeoutExpired: to = True
    wall = time.time() - t; row = dict(case=name, in_denominator=in_den, verdict='timeout' if to else 'no_output', passed=False, strict=False, engine_status='', objective='', reference_status='', reference_objective='', checker_objective='', rel_objective_error='', row_violation_rel='', bound_violation='', wall_s=round(wall, 3))
    st = json.load(open(js)) if os.path.exists(js) else {}
    if st: row['engine_status'] = st.get('status', ''); row['objective'] = st.get('objective', ''); row['verdict'] = st.get('status', 'no_output')
    A, c, rl, ru, l, u, const, names = mpsio.read_mps(path)
    rs, ro = highs(A, c, rl, ru, l, u, const); row['reference_status'] = rs; row['reference_objective'] = ro if ro is not None else ''
    if rs != 'Optimal': row['verdict'] = 'highs_failed'; return row
    if st.get('status') != 'optimal' or not os.path.exists(sol): return row
    x = np.zeros(len(names)); idx = {n: k for k, n in enumerate(names)}
    for line in open(sol):
        p = line.rstrip('\n').rsplit(None, 1)   # variable names may contain spaces (fixed-format MPS)
        if len(p) == 2 and p[0].strip() in idx: x[idx[p[0].strip()]] = float(p[1])
    fin = lambda z: np.where(np.isfinite(z), z, 0); Ax = A @ x
    rv = np.maximum(rl - Ax, 0) + np.maximum(Ax - ru, 0); rrel = float(np.max(rv / (1 + np.maximum(abs(fin(rl)), abs(fin(ru)))), initial=0))
    bv = float(np.max(np.maximum(l - x, 0) + np.maximum(x - u, 0), initial=0)); cobj = float(c @ x) + const
    tol = max(1e-6, 1e-7 * abs(ro)); e1 = abs(st['objective'] - ro); e2 = abs(cobj - st['objective']); rel = e1 / max(1.0, abs(ro))
    row.update(checker_objective=cobj, rel_objective_error=rel, row_violation_rel=rrel, bound_violation=bv)
    ok = e1 <= tol and e2 <= tol and rrel <= 1e-6 and bv <= 1e-6
    row['passed'] = bool(ok); row['strict'] = bool(ok and rrel <= 1e-8 and bv <= 1e-8 and rel <= 1e-8); row['verdict'] = 'pass' if ok else 'wrong_or_inaccurate'
    return row
for n in cases + extra:
    r = grade(n, n in cases); rows.append(r); print(n, r['verdict'], r['wall_s'], flush=True)
den = [r for r in rows if r['in_denominator']]; P = [r for r in den if r['passed'] and r['case'] not in fixed]; S = [r for r in P if r['strict']]
wrong = [r['case'] for r in den if r['verdict'] == 'wrong_or_inaccurate']; nonp = {}
for r in den:
    if not (r['passed'] and r['case'] not in fixed): nonp[r['case']] = r['verdict']
w = csv.DictWriter(open(os.path.join(a.out, 'ledger.csv'), 'w', newline=''), fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('RUN', a.label, '| cap', a.time_limit, 's | PASSES %d/%d' % (len(P), len(cases)), '| STRICT %d/%d' % (len(S), len(P)), '| WRONG', wrong, '| NOT PASSING', nonp)
for r in rows:
    if not r['in_denominator']: print('EXTRA', r['case'], r['verdict'], r['wall_s'], 'pass' if r['passed'] else 'no')
print('max row viol', max([r['row_violation_rel'] for r in P if r['row_violation_rel'] != ''] or [0]), 'max bound viol', max([r['bound_violation'] for r in P if r['bound_violation'] != ''] or [0]))
print('GATE_DONE', a.label)
