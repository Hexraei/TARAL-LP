"""Standalone MPS reader for the GPU measurement notebook. Fixed or free format, RANGES, bounds, objective RHS constant.
Returns A (csr), c, row lo/hi, col lo/hi, const, names."""
import numpy as np, scipy.sparse as sps
INF = np.inf
def _num(s): return float(s.replace('D', 'E').replace('d', 'e'))
def read_mps(path):
    lines = [l.rstrip('\n') for l in open(path)]
    fixed = False
    sec = None; rtype = {}; rorder = []; obj = None; cols = {}; corder = []; rhs = {}; rng = {}; bnd = {}
    def fields(line, sec, fixed):
        t = line.split()
        if not fixed: return t
        P = lambda a, b: line[a:b].strip()
        if sec == 'ROWS': return [P(1, 3), P(4, 12)]
        if sec == 'BOUNDS': return [x for x in (P(1, 3), P(4, 12), P(14, 22), P(24, 36)) if x != '']
        return [x for x in (P(4, 12), P(14, 22), P(24, 36), P(39, 47), P(49, 61)) if x != '']
    for attempt in (0, 1):
        fixed = bool(attempt)
        try:
            sec = None; rtype = {}; rorder = []; obj = None; cols = {}; corder = []; rhs = {}; rng = {}; bnd = {}
            for raw in lines:
                if not raw.strip() or raw[0] == '*': continue
                if raw[0] != ' ':
                    sec = raw.split()[0]; continue
                t = fields(raw, sec, fixed)
                if sec == 'ROWS':
                    if t[0] == 'N':
                        if obj is None: obj = t[1]
                    else: rtype[t[1]] = t[0]; rorder.append(t[1])
                elif sec == 'COLUMNS':
                    if "'MARKER'" in t: continue
                    j = t[0]
                    if j not in cols: cols[j] = {}; corder.append(j)
                    for k in range(1, len(t) - 1, 2): cols[j][t[k]] = cols[j].get(t[k], 0.0) + _num(t[k + 1])
                elif sec in ('RHS', 'RANGES'):
                    d = rhs if sec == 'RHS' else rng
                    rest = t[1:] if len(t) % 2 == 1 else t
                    for k in range(0, len(rest) - 1, 2): d[rest[k]] = _num(rest[k + 1])
                elif sec == 'BOUNDS':
                    typ = t[0]
                    if typ in ('FR', 'MI', 'PL', 'BV') and len(t) == 3: name, val = t[2], 0.0
                    elif len(t) == 4: name, val = t[2], _num(t[3])
                    elif len(t) == 3: name, val = t[1], _num(t[2])
                    else: name, val = t[-1], 0.0
                    lo, up = bnd.get(name, [0.0, INF])
                    if typ == 'LO': lo = val
                    elif typ == 'UP':
                        up = val
                        if val < 0 and lo == 0.0: lo = -INF
                    elif typ == 'FX': lo = up = val
                    elif typ == 'FR': lo, up = -INF, INF
                    elif typ == 'MI': lo = -INF
                    elif typ == 'PL': up = INF
                    elif typ == 'BV': lo, up = 0.0, 1.0
                    bnd[name] = [lo, up]
            break
        except (ValueError, IndexError):
            if attempt: raise
    rid = {r: i for i, r in enumerate(rorder)}; cid = {j: k for k, j in enumerate(corder)}
    m, n = len(rorder), len(corder); I, J, V = [], [], []; c = np.zeros(n)
    for j, col in cols.items():
        for r, a in col.items():
            if r == obj: c[cid[j]] += a
            elif r in rid: I.append(rid[r]); J.append(cid[j]); V.append(a)
    A = sps.csr_matrix((V, (I, J)), shape=(m, n)); lo = np.full(m, -INF); hi = np.full(m, INF)
    for r in rorder:
        i = rid[r]; k = rtype[r]; b = rhs.get(r, 0.0); R = rng.get(r)
        if k == 'L': hi[i] = b; lo[i] = b - abs(R) if R is not None else -INF
        elif k == 'G': lo[i] = b; hi[i] = b + abs(R) if R is not None else INF
        elif R is None: lo[i] = hi[i] = b
        elif R >= 0: lo[i], hi[i] = b, b + R
        else: lo[i], hi[i] = b + R, b
    l = np.array([bnd.get(j, [0.0, INF])[0] for j in corder]); u = np.array([bnd.get(j, [0.0, INF])[1] for j in corder])
    return A, c, lo, hi, l, u, -rhs.get(obj, 0.0), corder
