#!/usr/bin/env python3
"""Replay structural integer proofs from original MPS using independent HiGHS reader.
HiGHS is a reader only, never the proof oracle. Exact Fraction arithmetic.
Usage: structural_integer_replay.py MODEL.mps SOLVER.json
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from presolve_replay_check import load

def check(M,j):
    c=j['structural_certificate'];n=M['n'];m=M['m']
    if c['type']!='objective_lattice': assert all(M['isint']), 'noninteger column'
    rows=[{} for _ in range(m)]
    for k,col in enumerate(M['cols']):
        for i,v in col:rows[i][k]=F(v)
    def combination(indices):
        assert indices and len(indices)==len(set(indices))
        a=[F(0)]*n;b=F(0)
        for i in indices:
            assert isinstance(i,int) and 0<=i<m
            assert M['rlo'][i]==M['rup'][i] and math.isfinite(M['rlo'][i])
            rhs=F(M['rlo'][i]);assert rhs.denominator==1
            b+=rhs
            for k,v in rows[i].items():assert v.denominator==1;a[k]+=v
        return a,b
    if c['type']=='three_integer_congruence':
        assert M['n']==3 and M['m']==1 and not M['maximize'] and all(M['isint'])
        assert M['rlo'][0]==M['rup'][0]==0
        v=c['objective_col'];u,w=c['other_cols'];assert len({v,u,w})==3
        assert M['cost'][v]==1 and M['cost'][u]==M['cost'][w]==0
        assert M['lo'][v]==c['start']>=1 and M['lo'][u]==M['lo'][w]==0
        assert M['up'][u]==M['up'][w]==math.inf
        A=rows[0][v];B=-rows[0][u];C=-rows[0][w]
        assert all(a>0 and a.denominator==1 for a in (A,B,C))
        A,B,C=map(int,(A,B,C));inv=c['inverse'];assert B*inv%C==1
        end=c['end'];assert isinstance(end,int) and c['start']<=end<=1000000
        for z in range(c['start'],end):
            y=(A*z*inv)%C;assert B*y>A*z,'smaller feasible objective'
        x=list(map(F,j['x']));assert all(t.denominator==1 for t in x)
        assert x[v]==end and A*x[v]-B*x[u]-C*x[w]==0
        assert all(M['lo'][k]<=x[k]<=M['up'][k] for k in range(3))
        assert j['status']=='optimal' and j['has_solution']
        assert j['objective']==end+M['offset'] and j['best_bound']==j['objective'] and j['gap']==0
        return True
    if c['type']=='objective_lattice':
        t=c['objective_col'];q=c['denominator'];assert isinstance(q,int) and q>0
        cost=F(M['cost'][t]);sense=-1 if M['maximize'] else 1
        assert cost*sense<0 and not M['isint'][t] and M['lo'][t]==0
        assert all(v==0 for k,v in enumerate(M['cost']) if k!=t)
        found=set()
        for link in c['links']:
            r=link['row'];y=link['grid_col'];er=link['equality']
            assert r not in found;found.add(r)
            assert M['rlo'][r]==-math.inf and M['rup'][r]==0 and rows[r]=={t:F(1),y:F(-1)}
            assert not M['isint'][y] and M['lo'][y]==0 and M['up'][y]==math.inf
            assert M['rlo'][er]==M['rup'][er] and F(M['rlo'][er]).denominator==1
            assert rows[er][y]==q and all(M['isint'][k] and v.denominator==1 for k,v in rows[er].items() if k!=y)
            assert {i for i,v in M['cols'][y]}=={r,er}
        assert found=={i for i,v in M['cols'][t]} and found
        if math.isfinite(M['up'][t]):assert F(M['up'][t])*q==int(F(M['up'][t])*q)
        if j.get('has_solution'):
            x=list(map(F,j['x']));assert len(x)==n
            for k,v in enumerate(x):
                assert M['lo'][k]<=v<=M['up'][k]
                if M['isint'][k]:assert v.denominator==1
            for i in range(m):
                a=sum(v*x[k] for k,v in rows[i].items());assert M['rlo'][i]<=a<=M['rup'][i]
            obj=F(M['offset'])+sum(F(v)*x[k] for k,v in enumerate(M['cost']))
            assert abs(float(obj)-j['objective'])<1e-10
        return True
    if c['type']=='integer_parity_contradiction':
        assert j['status']=='infeasible' and not j.get('has_solution')
        a,b=combination(c['rows']);assert all(v%2==0 for v in a) and b%2==1
    elif c['type']=='integer_parity_unique':
        assert j['status']=='optimal' and j['has_solution']
        x=list(map(F,j['x']));assert len(x)==n
        known=set()
        for f in c['fixed_binary']:
            k=f['col'];assert k not in known and M['lo'][k]==0 and M['up'][k]==1
            assert f['value'] in (0,1)
            a,b=combination(f['rows']);assert a[k]%2==1 and all(a[t]%2==0 for t in range(n) if t!=k)
            assert b%2==f['value'] and x[k]==f['value'];known.add(k)
        assert known=={k for k in range(n) if M['lo'][k]==0 and M['up'][k]==1}
        # Every remaining column must have one equality incident to it, with all
        # other columns already unique. Thus the full integer point is unique.
        for k in range(n):
            if k in known:continue
            assert len(M['cols'][k])==1
            i,v=M['cols'][k][0]
            assert v!=0 and M['rlo'][i]==M['rup'][i]
            assert all(t in known for t in rows[i] if t!=k)
        for k,v in enumerate(x):
            assert v.denominator==1 and M['lo'][k]<=v<=M['up'][k]
        for i in range(m):
            act=sum(v*x[k] for k,v in rows[i].items())
            assert M['rlo'][i]<=act<=M['rup'][i]
        obj=F(M['offset'])+sum(F(v)*x[k] for k,v in enumerate(M['cost']))
        assert abs(float(obj)-j['objective'])<=1e-12*max(1,abs(float(obj)))
        assert j['best_bound']==j['objective'] and j['gap']==0
    else:raise AssertionError('unsupported proof type')
    return True
if __name__=='__main__':
    check(load(sys.argv[1]),json.load(open(sys.argv[2])));print('PASS exact structural proof replay')
