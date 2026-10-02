"""Mittelmann-set LP runs: src/ dual (headline) and primal simplex vs HiGHS reference on the ORIGINAL MPS. Public inputs, hashes recorded. Measured."""
import sys, os, json, time, subprocess, hashlib, bz2, gzip, shutil, urllib.request
import numpy as np, scipy.sparse as sps
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import mpsio
from scipy.optimize import linprog
KIT, GROUP = sys.argv[1], sys.argv[2]
TLC, TLH = float(sys.argv[3]), float(sys.argv[4])
def run(cmd, tl, tag):
    sol = '/tmp/%s.sol' % tag; js = '/tmp/%s.json' % tag
    for f in (sol, js):
        if os.path.exists(f): os.remove(f)
    t = time.time()
    try: subprocess.run(cmd + ['--sol', sol, '--json', js], timeout=tl + 30, capture_output=True)
    except subprocess.TimeoutExpired: pass
    return time.time() - t, (json.load(open(js)) if os.path.exists(js) else {}), sol
def check(M, sol, ref):
    A, c, rl, ru, l, u, const, names = M; idx = {n: k for k, n in enumerate(names)}; x = np.zeros(len(names))
    if not os.path.exists(sol): return None
    for line in open(sol):
        p = line.rstrip('\n').rsplit(None, 1)
        if len(p) == 2 and p[0].strip() in idx: x[idx[p[0].strip()]] = float(p[1])
    fin = lambda z: np.where(np.isfinite(z), z, 0); Ax = A @ x
    rv = np.maximum(rl - Ax, 0) + np.maximum(Ax - ru, 0); rr = float(np.max(rv / (1 + np.maximum(abs(fin(rl)), abs(fin(ru)))), initial=0))
    bv = float(np.max(np.maximum(l - x, 0) + np.maximum(x - u, 0), initial=0)); obj = float(c @ x) + const
    err = abs(obj - ref) / (1 + abs(ref)); return dict(objective=obj, rel_row=float('%.3g' % rr), bound=float('%.3g' % bv), obj_rel_err=float('%.3g' % err), meets_1e6=bool(rr <= 1e-6 and bv <= 1e-6 and err <= 1e-6))
def highs(M):
    A, c, rl, ru, l, u, const, names = M; fl = np.isfinite(rl); fh = np.isfinite(ru); eq = fl & fh & (rl == ru); ih = fh & ~eq; il = fl & ~eq
    Aub = sps.vstack([A[ih], -A[il]]).tocsr(); bub = np.concatenate([ru[ih], -rl[il]]); t = time.time()
    r = linprog(c, A_ub=Aub if Aub.shape[0] else None, b_ub=bub if Aub.shape[0] else None, A_eq=A[eq] if eq.sum() else None, b_eq=rl[eq] if eq.sum() else None,
                bounds=list(zip(np.where(np.isfinite(l), l, None), np.where(np.isfinite(u), u, None))), method='highs', options={'time_limit': TLH})
    return (float(r.fun) + const if r.status == 0 else None), time.time() - t
SRC = {
 'fome11':('https://plato.asu.edu/ftp/lptestset/fome/fome11.bz2','bz2'), 'fome12':('https://plato.asu.edu/ftp/lptestset/fome/fome12.bz2','bz2'),
 'fome13':('https://plato.asu.edu/ftp/lptestset/fome/fome13.bz2','bz2'), 'fome21':('https://plato.asu.edu/ftp/lptestset/fome/fome21.bz2','bz2'),
 'rail507':('https://plato.asu.edu/ftp/lptestset/rail/rail507.bz2','bz2'), 'rail516':('https://plato.asu.edu/ftp/lptestset/rail/rail516.bz2','bz2'),
 'rail582':('https://plato.asu.edu/ftp/lptestset/rail/rail582.bz2','bz2'), 'rail2586':('https://plato.asu.edu/ftp/lptestset/rail/rail2586.bz2','bz2'),
 'rail4284':('https://plato.asu.edu/ftp/lptestset/rail/rail4284.bz2','bz2'),
}
for n in ('pds-20','pds-30','pds-40','pds-50','pds-60','pds-70','pds-80','pds-90','pds-100'): SRC[n] = ('https://plato.asu.edu/ftp/lptestset/pds/%s.bz2' % n, 'bz2')
for n in ('cre-a','cre-b','cre-c','cre-d','ken-07','ken-11','ken-13','ken-18','osa-07','osa-14','osa-30','osa-60','pds-02','pds-06','pds-10'): SRC[n] = ('https://www.netlib.org/lp/data/kennington/%s.gz' % n, 'gz')
GROUPS = {'H2A': ['rail2586','rail4284'],'N1': ['fome11', 'fome12'],'N2': ['fome13', 'fome21'],'N3': ['rail2586'],'N4': ['rail4284'],'N5': ['rail507', 'rail516', 'rail582', 'ken-18'],'N6': ['pds-20', 'pds-30', 'pds-40', 'pds-50'],'N7': ['pds-60', 'pds-70', 'pds-80', 'pds-90', 'pds-100'],'N8': ['ken-07', 'cre-a', 'cre-b', 'osa-07', 'pds-02', 'ken-11', 'cre-c'],'N9': ['cre-d', 'osa-14', 'pds-06', 'ken-13', 'osa-30', 'pds-10', 'osa-60'],'A1': ['fome11','fome12','fome13','fome21','rail507','rail516','rail582','rail2586','rail4284'],
          'A2': ['ken-07','cre-a','cre-b','osa-07','pds-02','ken-11','cre-c','cre-d','osa-14','pds-06','ken-13','osa-30','pds-10','osa-60','ken-18'],
          'R1': ['rail507','rail516','rail582'], 'R2': ['rail2586'], 'R3': ['rail4284'], 'F': ['fome21'], 'P1': ['pds-100','pds-90'], 'P2': ['pds-80','pds-70'],
          'B': ['pds-20','pds-30','pds-40','pds-50','pds-60','pds-70','pds-80','pds-90','pds-100']}
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
MITDIR=os.environ.get('MITDIR','/app/mitwork'); os.makedirs(MITDIR, exist_ok=True)
for name in GROUPS[GROUP]:
    url, kind = SRC[name]; d = MITDIR + '/' + name; os.makedirs(d, exist_ok=True); rec = dict(case=name, url=url)
    try:
        raw = d + '/' + name + '.' + kind; open(raw, 'wb').write(urllib.request.urlopen(url, timeout=300).read()); rec['download_sha256'] = sha(raw); rec['download_bytes'] = os.path.getsize(raw)
        mpc = d + '/' + name + '.mpc'; open(mpc, 'wb').write((bz2 if kind == 'bz2' else gzip).open(raw).read())
        mps = d + '/' + name + '.mps'; subprocess.run('%s/emps %s > %s' % (KIT, mpc, mps), shell=True, check=True)
        rec['mps_sha256'] = sha(mps); rec['mps_bytes'] = os.path.getsize(mps)
        M = mpsio.read_mps(mps); A = M[0]; rec.update(rows=int(A.shape[0]), cols=int(A.shape[1]), nnz=int(A.nnz))
    except Exception as e:
        rec['error'] = 'input: ' + repr(e)[:200]; print('MIT_ROW', json.dumps(rec), flush=True); continue
    try: ref, th = (None, 0.0) if TLH <= 0 else highs(M)
    except Exception as e: ref, th = None, 0.0
    rec.update(reference=ref, highs_wall=round(th, 2)); print('CASE', name, rec['rows'], rec['cols'], rec['nnz'], 'HiGHS ref', ref, 'wall', round(th, 1), flush=True)
    for meth in ('dual', 'primal'):
        cmd = [KIT + '/taral', mps, '--time-limit', str(TLC)] + (['--method', 'dual'] if meth == 'dual' else [])
        sol = d + '/' + meth + '.sol'; js = d + '/' + meth + '.json'; t = time.time()
        try: subprocess.run(cmd + ['--sol', sol, '--json', js], timeout=TLC + 60, capture_output=True)
        except subprocess.TimeoutExpired: pass
        w = time.time() - t; st = json.load(open(js)) if os.path.exists(js) else {}
        ck = check(M, sol, ref) if ref is not None else None
        obj = st.get('objective'); ok = ref is not None and st.get('status') == 'optimal' and obj is not None and abs(obj - ref) <= 1e-6 * (1 + abs(ref))
        r = dict(wall=round(w, 2), status=st.get('status', 'no_output'), iterations=st.get('iterations'), objective=obj, obj_matches_ref=bool(ok)); r.update(ck or {})
        rec[meth] = r; print('  ', meth.upper(), r, flush=True)
    print('MIT_ROW', json.dumps(rec), flush=True)
print('MIT_DONE', GROUP)
