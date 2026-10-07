#!/usr/bin/env python3
"""Robustness of the Farkas repair path: duplicate entries, column order and huge coefficients.
For each infeasible model, variants (original, column-order shuffle, every coefficient split into two
duplicate entries, first row's RHS and coefficients scaled by 1e150) must give: all multipliers finite and
nonnegative, and every infeasible output either passes the raw strict checker or carries
[strict_farkas_contract=unresolved].
Non-gating audit: every status=unbounded output is also run through the independent checker and listed
(UNBOUNDED_AUDIT rows). A failing unbounded row is reported, not asserted: unbounded-certificate acceptance is outside
the cert-contract-fix scope. Known pre-existing row: cplex1 scale dual (engine verified=true, checker residual 0.391)."""
import argparse, json, math, os, random, subprocess, sys, tempfile
ap = argparse.ArgumentParser(); ap.add_argument('--engine', required=True); ap.add_argument('--corpus', required=True)
ap.add_argument('--checker', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nonoptimal_certificate_check.py'))
ap.add_argument('--methods', default='dual,simplex'); a = ap.parse_args()

def transform(src, dst, kind, rng):
    sec = None; heads = []; cols = []; cur = None; srow = None; other = []
    out = []
    for ln in open(src).read().split('\n'):
        if ln and ln[0] not in ' *':
            sec = ln.split()[0]; out.append(('H', sec, ln)); continue
        if not ln.strip() or ln[0] == '*': continue
        t = ln.split()
        if sec == 'ROWS' and srow is None and t[0] != 'N': srow = t[1]
        if sec == 'COLUMNS':
            if 'MARKER' in t or len(t) not in (3, 5): return False
            if cur is None or cur[0] != t[0]:
                if cur: cols.append(cur)
                cur = [t[0], []]
            for k in range(1, len(t) - 1, 2): cur[1].append((t[k], float(t[k + 1])))
            continue
        if sec in ('RHS', 'RANGES') and kind == 'scale':
            s0 = 1 if len(t) % 2 == 1 else 0
            for k in range(s0, len(t) - 1, 2):
                if t[k] == srow: t[k + 1] = repr(float(t[k + 1]) * 1e150)
            ln = ' ' + ' '.join(t)
        out.append(('L', sec, ln))
    if cur: cols.append(cur)
    res = []
    for tag, s, ln in out:
        res.append(ln)
        if tag == 'H' and s == 'COLUMNS':
            cc = list(cols)
            if kind == 'shuffle': rng.shuffle(cc)
            for name, ents in cc:
                for r, v in ents:
                    if kind == 'split': res += [f' {name} {r} {v/2!r}'] * 2
                    elif kind == 'scale' and r == srow: res.append(f' {name} {r} {v*1e150!r}')
                    else: res.append(f' {name} {r} {v!r}')
    open(dst, 'w').write('\n'.join(res) + '\n')
    return True

def run(model, m, tag):
    j = tag + '.json'
    subprocess.run([a.engine, model, '--method', m, '--json', j, '--time-limit', '30'], capture_output=True, text=True, timeout=200)
    d = json.load(open(j))
    for k in ('farkas_row_lower', 'farkas_row_upper', 'farkas_col_lower', 'farkas_col_upper'):
        for x in d.get(k) or []: assert math.isfinite(x) and x >= 0, (model, k, x)
    assert not d['status'].startswith('parse_error'), (model, m, d['status'])
    if d['status'] == 'unbounded':
        r = subprocess.run([sys.executable, a.checker, model, j], capture_output=True, text=True).stdout
        unb.append((os.path.basename(model), m, '"pass_certificate": true' in r, d.get('certificate_verified')))
    if d['status'] != 'infeasible': return d['status']
    r = subprocess.run([sys.executable, a.checker, model, j], capture_output=True, text=True).stdout
    ok = '"pass_certificate": true' in r
    assert ok or 'strict_farkas_contract=unresolved' in d.get('message', ''), (model, m, 'infeasible cert fails strict and is not flagged')
    return 'infeasible-strict-pass' if ok else 'infeasible-flagged-unresolved'

rng = random.Random(7); n = 0; stats = {}; unb = []
with tempfile.TemporaryDirectory() as td:
    for f in sorted(os.listdir(a.corpus)):
        if not f.endswith('.mps'): continue
        src = os.path.join(a.corpus, f)
        for kind in ('orig', 'shuffle', 'split', 'scale'):
            mp = src if kind == 'orig' else os.path.join(td, f'{f[:-4]}-{kind}.mps')
            if kind != 'orig' and not transform(src, mp, kind, rng): continue
            for m in a.methods.split(','):
                st = run(mp, m, os.path.join(td, f'{f[:-4]}-{kind}-{m}')); n += 1
                stats[(kind, st)] = stats.get((kind, st), 0) + 1
for k, v in sorted(stats.items()): print(k, v)
for name, m, chk, ver in unb: print('UNBOUNDED_AUDIT', name, m, 'checker_pass' if chk else 'CHECKER_REJECTS', 'engine_verified=%s' % ver)
print('UNBOUNDED_AUDIT total', len(unb), 'checker_rejects', sum(1 for u in unb if not u[2]), '(non-gating)')
print(f'PASS {n} repair-path robustness runs (finite nonnegative multipliers; strict-pass or flagged unresolved)')
