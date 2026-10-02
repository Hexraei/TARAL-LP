#!/usr/bin/env python3
"""Independent original-MPS LP Farkas/recession proof checker.
No solver parser, presolved model, reported residual, or verified flag is trusted.
Farkas inequalities: -A'x row_lower + A'x row_upper - col_lower + col_upper = 0;
weighted bound contradiction strictly positive. All multipliers nonnegative, L1=1.
Unbounded: feasible original anchor plus normalized recession ray and improving slope.
Finite precision certificates, not exact rational proofs. Stationarity <=1e-12,
primal/recession residual <=1e-8, normalized strict margin >1e-8.
Integer and quadratic models are rejected. Original first-appearance array order.
"""
import argparse
import json
import math
from pathlib import Path

INF = float('inf')


class UnsupportedModel(ValueError):
    """Not an LP; do not retry these files as fixed-column MPS."""


def read_model(path):
    def parse(fixed):
        section = None
        rows, columns, costs, bounds, rhs, ranges, sets = {}, {}, {}, {}, {}, {}, {}
        objective, offset, sense = None, 0.0, 1
        def number(value):
            v = float(value.replace('D', 'E').replace('d', 'e'))
            if not math.isfinite(v):
                raise ValueError('nonfinite MPS coefficient')
            return v
        for raw in Path(path).read_text().splitlines():
            if not raw.strip() or raw.startswith('*'):
                continue
            t = raw.split()
            if not raw[0].isspace():
                section = t[0]
                if section == 'ENDATA':
                    break
                if section not in ('NAME', 'OBJSENSE', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS'):
                    raise UnsupportedModel('unsupported section: ' + section)
                if section != 'OBJSENSE' or len(t) < 2:
                    continue
                t = t[1:]
            elif fixed and section in ('ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS'):
                f = lambda a, b: raw[a:b].strip()
                if section == 'ROWS':
                    t = [f(1, 3), f(4, 12)]
                elif section == 'BOUNDS':
                    t = [f(1, 3), f(4, 12), f(14, 22)] + ([f(24, 36)] if f(24, 36) else [])
                else:
                    t = [f(4, 12), f(14, 22), f(24, 36)] + ([f(39, 47), f(49, 61)] if f(39, 47) else [])
            if section == 'OBJSENSE':
                if t[0] not in ('MIN', 'MINIMIZE', 'MAX', 'MAXIMIZE'):
                    raise ValueError('bad objective sense')
                sense = -1 if t[0].startswith('MAX') else 1
            elif section == 'ROWS':
                if len(t) != 2 or t[0] not in ('N', 'E', 'L', 'G'):
                    raise ValueError('bad row')
                if t[0] == 'N' and objective is None:
                    objective = t[1]
                rows[t[1]] = t[0]
            elif section == 'COLUMNS':
                if "'MARKER'" in t:
                    raise UnsupportedModel('integer markers are not LP certificates')
                if len(t) not in (3, 5):
                    raise ValueError('bad column')
                j = t[0]
                col = columns.setdefault(j, {})
                costs.setdefault(j, 0.0)
                for k in range(1, len(t), 2):
                    r, v = t[k], number(t[k + 1])
                    if r not in rows:
                        raise ValueError('unknown row')
                    if r == objective:
                        costs[j] += v
                    elif rows[r] != 'N':
                        col[r] = col.get(r, 0.0) + v
            elif section in ('RHS', 'RANGES'):
                start = 1 if fixed or len(t) % 2 else 0
                if start and t[0]:
                    if sets.setdefault(section, t[0]) != t[0]:
                        continue
                for k in range(start, len(t), 2):
                    r, v = t[k], number(t[k + 1])
                    if r not in rows:
                        raise ValueError('unknown rhs/range row')
                    if section == 'RHS' and r == objective:
                        offset = -v
                    else:
                        (rhs if section == 'RHS' else ranges)[r] = v
            elif section == 'BOUNDS':
                typ = t[0]
                if typ in ('BV', 'LI', 'UI'):
                    raise UnsupportedModel('integer bounds are not LP certificates')
                no_value = typ in ('FR', 'MI', 'PL')
                if fixed or len(t) == 4 or (len(t) == 3 and no_value):
                    name, j = t[1:3]
                    value = t[3] if len(t) == 4 else None
                else:
                    name, j = '', t[1]
                    value = t[2] if len(t) == 3 else None
                if name and sets.setdefault('BOUNDS', name) != name:
                    continue
                if j not in columns:
                    raise ValueError('unknown bound column')
                v = number(value) if value is not None else 0.0
                lo, up = bounds.get(j, (0.0, INF))
                if typ == 'LO': lo = v
                elif typ == 'UP':
                    up = v
                    if v < 0 and lo == 0: lo = -INF
                elif typ == 'FX': lo = up = v
                elif typ == 'FR': lo, up = -INF, INF
                elif typ == 'MI': lo = -INF
                elif typ == 'PL': up = INF
                else: raise ValueError('unsupported bound: ' + typ)
                bounds[j] = (lo, up)
        row_bounds = {}
        for r, kind in rows.items():
            if kind == 'N': continue
            b, width = rhs.get(r, 0.0), ranges.get(r)
            lo, up = (b, b) if kind == 'E' else ((-INF, b) if kind == 'L' else (b, INF))
            if width is not None:
                if kind == 'L': lo = b - abs(width)
                elif kind == 'G': up = b + abs(width)
                elif width >= 0: up = b + width
                else: lo = b + width
            row_bounds[r] = (lo, up)
        return columns, costs, bounds, row_bounds, offset, sense
    try:
        return parse(False)
    except UnsupportedModel:
        raise
    except (ValueError, IndexError):
        return parse(True)


def verify(path, certificate):
    cols, costs, bounds, rows, offset, sense = read_model(path)
    def vector(key, size, positive=False):
        v = certificate[key]
        if not isinstance(v,list) or len(v)!=size or any(not isinstance(a,(float,int)) or isinstance(a,bool) or not math.isfinite(a) or (positive and a<0) for a in v):
            raise ValueError('invalid '+key)
        return v
    n,m=len(cols),len(rows)
    residual=margin=0.0
    if certificate.get('status')=='infeasible':
        rl=vector('farkas_row_lower',m,True);ru=vector('farkas_row_upper',m,True)
        cl=vector('farkas_col_lower',n,True);cu=vector('farkas_col_upper',n,True)
        norm=math.fsum(rl+ru+cl+cu); terms=[]
        def term(a,b,sgn):
            if a:
                if not math.isfinite(b): raise ValueError('multiplier on infinite bound')
                terms.append(sgn*a*b)
        for (r,(l,u)),a,b in zip(rows.items(),rl,ru):term(a,l,1);term(b,u,-1)
        y=dict(zip(rows,[b-a for a,b in zip(rl,ru)]))
        for (j,col),a,b in zip(cols.items(),cl,cu):
            z=math.fsum([b,-a]+[v*y[r] for r,v in col.items()])
            residual=max(residual,abs(z))
            l,u=bounds.get(j,(0.,INF));term(a,l,1);term(b,u,-1)
            if z:
                bound=l if z>0 else u
                if not math.isfinite(bound):
                    for r,coeff in col.items():
                        if coeff and all(k==j or other.get(r,0.)==0 for k,other in cols.items()):
                            endpoint=rows[r][0 if (z>0)==(coeff>0) else 1]
                            candidate=endpoint/coeff
                            if math.isfinite(candidate):
                                bound=candidate if not math.isfinite(bound) else (max(bound,candidate) if z>0 else min(bound,candidate))
                if not math.isfinite(bound): raise ValueError('uncancelled stationarity on unbounded column')
                terms.append(z*bound)
        margin=math.fsum(terms)/(1+math.fsum(abs(t) for t in terms))
        ok=abs(norm-1)<=1e-8 and residual<=1e-12 and margin>1e-8
    elif certificate.get('status')=='unbounded':
        x=vector('x',n);d=vector('ray',n);act={r:[] for r in rows};direction={r:[] for r in rows}
        def check(v,w,l,u):
            nonlocal residual
            if math.isfinite(l):residual=max(residual,max(0.,(l-v)/(1+abs(l))),max(0.,-w))
            if math.isfinite(u):residual=max(residual,max(0.,(v-u)/(1+abs(u))),max(0.,w))
        for (j,col),v,w in zip(cols.items(),x,d):
            check(v,w,*bounds.get(j,(0.,INF)))
            for r,a in col.items():act[r].append(a*v);direction[r].append(a*w)
        for r,(l,u) in rows.items():check(math.fsum(act[r]),math.fsum(direction[r]),l,u)
        margin=-sense*math.fsum(costs[j]*w for j,w in zip(cols,d))/(1+math.fsum(abs(v) for v in costs.values()))
        ok=residual<=1e-8 and abs(max(map(abs,d),default=0.)-1)<=1e-8 and margin>1e-8
    else: ok=False
    return dict(pass_certificate=bool(ok and math.isfinite(residual) and math.isfinite(margin)),residual=residual,margin=margin)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('model');ap.add_argument('certificate');a=ap.parse_args()
    try:r=verify(a.model,json.loads(Path(a.certificate).read_text()))
    except (ValueError,TypeError,KeyError,IndexError,OSError,OverflowError) as e:r=dict(pass_certificate=False,error=str(e))
    print(json.dumps(r,allow_nan=False));raise SystemExit(0 if r['pass_certificate'] else 1)
if __name__=='__main__':main()
