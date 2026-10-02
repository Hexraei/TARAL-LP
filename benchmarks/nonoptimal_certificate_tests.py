#!/usr/bin/env python3
"""Synthetic original-model proofs, SciPy HiGHS reference and corrupt-ray rejection."""
import argparse, copy, json, random, subprocess, tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from nonoptimal_certificate_check import verify,read_model


def mps(A,lo,up,c,bounds,maximize=False):
    lines=['NAME PROOF','OBJSENSE '+('MAX' if maximize else 'MIN'),'ROWS',' N obj']
    for i in range(len(A)):
        if np.isfinite(lo[i]):lines.append(f' G l{i}')
        if np.isfinite(up[i]):lines.append(f' L u{i}')
    lines+=['COLUMNS']
    for j in range(len(c)):
        lines.append(f' x{j} obj {c[j]:.17g}')
        for i in range(len(A)):
            for prefix,finite in [('l',np.isfinite(lo[i])),('u',np.isfinite(up[i]))]:
                if finite and A[i][j]:lines.append(f' x{j} {prefix}{i} {A[i][j]:.17g}')
    lines+=['RHS']
    for i in range(len(A)):
        if np.isfinite(lo[i]):lines.append(f' rhs l{i} {lo[i]:.17g}')
        if np.isfinite(up[i]):lines.append(f' rhs u{i} {up[i]:.17g}')
    lines+=['BOUNDS']
    for j,(l,u) in enumerate(bounds):
        lines.append(f' FR b x{j}')
        if np.isfinite(l):lines.append(f' LO b x{j} {l:.17g}')
        if np.isfinite(u):lines.append(f' UP b x{j} {u:.17g}')
    return '\n'.join(lines+['ENDATA',''])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--engine',required=True);a=ap.parse_args()
    total=rejected=0;rng=random.Random(57463)
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/'model.mps';js=Path(td)/'result.json'
        def run(text,status):
            nonlocal total,rejected
            path.write_text(text)
            cols,costs,bnd,rows,_,sense=read_model(path)
            mat=[];rhs=[]
            for r,(l,u) in rows.items():
                row=[cols[j].get(r,0.) for j in cols]
                if np.isfinite(l):mat.append([-v for v in row]);rhs.append(-l)
                if np.isfinite(u):mat.append(row);rhs.append(u)
            reference=linprog([sense*costs[j] for j in cols],A_ub=mat or None,b_ub=rhs or None,
                             bounds=[bnd.get(j,(0.,np.inf)) for j in cols],method='highs',options={'presolve':False})
            assert reference.status=={'infeasible':2,'unbounded':3}[status],reference.message
            for method in ['simplex','dual']:
                subprocess.run([a.engine,str(path),'--method',method,'--json',str(js)],check=True,capture_output=True,timeout=15)
                cert=json.loads(js.read_text());assert cert['status']==status,(text,method,cert)
                assert cert['certificate_verified'] and verify(path,cert)['pass_certificate'],cert
                total+=1
                keys=['farkas_row_lower','farkas_row_upper','farkas_col_lower','farkas_col_upper'] if status=='infeasible' else ['x','ray']
                for key in keys:
                    if not cert[key]:continue
                    for v in [None,float('nan'),float('inf'),True]:
                        bad=copy.deepcopy(cert);bad[key][0]=v
                        try:passed=verify(path,bad)['pass_certificate']
                        except (ValueError,TypeError):passed=False
                        assert not passed;rejected+=1
                    bad=copy.deepcopy(cert);bad[key]=bad[key][:-1]
                    try:passed=verify(path,bad)['pass_certificate']
                    except ValueError:passed=False
                    assert not passed;rejected+=1
                bad=copy.deepcopy(cert)
                if status=='infeasible':
                    for k in keys:bad[k]=[0.]*len(bad[k])
                else:bad['ray']=[-v for v in bad['ray']]
                assert not verify(path,bad)['pass_certificate'];rejected+=1
                # Mutating original RHS/bounds must invalidate a formerly valid proof.
                saved=path.read_text()
                if status=='infeasible':
                    altered=mps([],[],[],[0.]*len(cols),[(-np.inf,np.inf)]*len(cols))
                else:altered=mps([],[],[],[0.]*len(cols),[(0.,0.)]*len(cols))
                path.write_text(altered)
                try:passed=verify(path,cert)['pass_certificate']
                except ValueError:passed=False
                assert not passed;rejected+=1;path.write_text(saved)
        # Independently exercise compensated near-zero stationarity with a
        # free direct bound and original singleton-row bound evidence.
        path.write_text(mps([[1.],[1.]],[2.,-np.inf],[np.inf,1.],[0.],[(-np.inf,np.inf)]))
        proof=dict(status='infeasible',farkas_row_lower=[0.5+1e-14,0.],
                   farkas_row_upper=[0.,0.5-1e-14],farkas_col_lower=[0.],farkas_col_upper=[0.])
        assert verify(path,proof)['pass_certificate'];total+=1
        path.write_text(mps([[1.,1.],[1.,1.]],[2.,-np.inf],[np.inf,1.],[0.,0.],[(-np.inf,np.inf)]*2))
        proof['farkas_col_lower']=[0.,0.];proof['farkas_col_upper']=[0.,0.]
        try:passed=verify(path,proof)['pass_certificate']
        except ValueError:passed=False
        assert not passed;rejected+=1
        run(mps([],[],[],[0.],[(2.,1.)]),'infeasible')
        run(mps([[1.]],[2.],[1.],[0.],[(-np.inf,np.inf)]),'infeasible')
        run(mps([[1.,1.]],[3.],[np.inf],[1.,0.],[(0.,1.),(0.,1.)]),'infeasible')
        run(mps([],[],[],[-1.],[(0.,np.inf)]),'unbounded')
        run(mps([],[],[],[1.],[(-np.inf,np.inf)],True),'unbounded')
        run(mps([],[],[],[1.],[(0.,np.inf)],True),'unbounded')
        run(mps([],[],[],[-1.],[(-np.inf,0.)],True),'unbounded')
        run(mps([],[],[],[1.],[(2.,1.)],True),'infeasible')
        run(mps([[1.,-1.]],[0.],[0.],[-1.,0.],[(0.,np.inf)]*2),'unbounded')
        for t in range(60):
            n=rng.randint(1,6);v=[rng.randint(1,5) for _ in range(n)]
            run(mps([v,v],[1.,-np.inf],[np.inf,0.],[rng.randint(-3,3) for _ in range(n)],[(-np.inf,np.inf)]*n,t%2==1),'infeasible')
            # Coupled recession ray with a feasible nonzero equality anchor.
            A=[[1.,-1.]+[0.]*(n-1)];c=[-1.,0.]+[0.]*(n-1)
            if t%2:c=[-v for v in c]
            run(mps(A,[float(t%5)],[float(t%5)],c,[(0.,np.inf)]*(n+1),t%2==1),'unbounded')
    print(f'PASS {total} original-model proof/reference checks; {rejected} corrupt/nonfinite/model-mutation rejections')
if __name__=='__main__':main()
