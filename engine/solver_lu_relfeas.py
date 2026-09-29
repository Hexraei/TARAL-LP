"""Experimental from-scratch two-phase revised simplex, for small standard-form LPs.
No optimization solver is used. NumPy handles matrix algebra only.
"""
from pathlib import Path
import time
import numpy as np
from scipy.linalg import lu_factor,lu_solve


def read_mps(path):
    section = None; rows = {}; columns = {}; rhs = {}; objective = None
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('*'): continue
        tok = line.split()
        if tok[0] in {'NAME','ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA'} and (line[0] != ' '):
            section = tok[0]
            if section in {'RANGES','BOUNDS'}: raise NotImplementedError(section)
            continue
        if section == 'ROWS':
            kind, name = tok[:2]
            if kind == 'N':
                if objective is None: objective = name
            else: rows[name] = kind
        elif section == 'COLUMNS':
            name = tok[0]
            if name not in columns: columns[name] = {}
            for k in range(1,len(tok),2):
                columns[name][tok[k]] = columns[name].get(tok[k],0.)+float(tok[k+1].replace('D','E'))
        elif section == 'RHS':
            for k in range(1 if len(tok)%2 else 0,len(tok),2): rhs[tok[k]] = float(tok[k+1].replace('D','E'))
    names = list(columns); rnames = list(rows)
    A=np.array([[columns[x].get(r,0.) for x in names] for r in rnames],dtype=float)
    b=np.array([rhs.get(r,0.) for r in rnames],dtype=float)
    c=np.array([columns[x].get(objective,0.) for x in names],dtype=float)
    return A,b,c,list(rows.values()),names


def _revised(A,b,c,basis,allowed=None,max_iter=10000,tol=1e-9):
    m,n=A.shape; basis=list(basis); allowed=set(range(n)) if allowed is None else set(allowed)
    for it in range(max_iter):
        B=A[:,basis]
        lu,piv=lu_factor(B,check_finite=False)
        xB=lu_solve((lu,piv),b,check_finite=False)
        if np.min(xB)<-max(1e-7,1e-12*np.max(np.abs(b))): raise ArithmeticError('lost primal feasibility')
        pi=lu_solve((lu,piv),c[basis],trans=1,check_finite=False); reduced=c-A.T@pi
        candidate=min((j for j in sorted(allowed) if j not in basis and reduced[j]<-tol),key=lambda j:(reduced[j],j),default=None)
        if candidate is None:
            x=np.zeros(n); x[basis]=xB
            return x,float(c@x),basis,it
        d=lu_solve((lu,piv),A[:,candidate],check_finite=False)
        pivot_cutoff=max(1e-8,1e-7*np.max(np.abs(d)))
        ratios=[(max(0.,xB[i])/d[i], basis[i],i) for i in range(m) if d[i]>pivot_cutoff]
        if not ratios: raise ArithmeticError('unbounded')
        _,_,leave=min(ratios)
        basis[leave]=candidate
    raise ArithmeticError('iteration limit')


def solve_lp(A,b,c,kinds):
    """min c*x with rows <=, >= or = and x>=0. Phase I artificials, phase II original."""
    A=np.asarray(A,float).copy(); b=np.asarray(b,float).copy(); c=np.asarray(c,float)
    kinds=list(kinds); m,n=A.shape
    for i in range(m):
        if b[i]<0:
            A[i]*=-1; b[i]*=-1
            kinds[i]={'L':'G','G':'L','E':'E'}[kinds[i]]
    mats=[A]; basis=[]; artificial=[]; extras=[]
    for i,kind in enumerate(kinds):
        if kind not in {'L','G','E'}: raise ValueError(kind)
        if kind!='E':
            v=np.zeros(m); v[i]=1 if kind=='L' else -1
            extras.append(v); slack=n+len(extras)-1
        if kind=='L': basis.append(slack)
        else: basis.append(None)
    if extras: mats.append(np.column_stack(extras))
    core=np.column_stack(mats); astart=core.shape[1]
    for i,j in enumerate(basis):
        if j is None:
            v=np.zeros(m); v[i]=1
            artificial.append(v); basis[i]=astart+len(artificial)-1
    full=np.column_stack([core]+([np.column_stack(artificial)] if artificial else []))
    phase1=np.zeros(full.shape[1]); phase1[astart:]=1.
    x,penalty,basis,it1=_revised(full,b,phase1,basis)
    if penalty>1e-7: raise ArithmeticError(f'infeasible; phase-I penalty={penalty}')
    # Pivot zero-valued artificial basic variables out. If a row is dependent,
    # drop it for phase II (but retain original rows for residual checks).
    keep=list(range(m))
    for i in range(m-1,-1,-1):
        if basis[i] < astart: continue
        B=full[np.ix_(keep,basis)]
        for j in range(astart):
            if j in basis: continue
            d=np.linalg.solve(B,full[np.ix_(keep,[j])].ravel())
            if abs(d[i])>1e-8:
                basis[i]=j; break
        if basis[i]>=astart:
            keep.pop(i); basis.pop(i)
    phase2=np.r_[c,np.zeros(astart-n)]
    x,z,basis,it2=_revised(core[keep],b[keep],phase2,basis)
    return {'x':x[:n], 'objective':z, 'phase1_iterations':it1,'phase2_iterations':it2,
            'kept_rows':len(keep),'max_residual':float(max([0.]+[
                 (max(0.,A[i]@x[:n]-b[i]) if kinds[i]=='L' else
                  max(0.,b[i]-A[i]@x[:n]) if kinds[i]=='G' else abs(A[i]@x[:n]-b[i]))
                 for i in range(m)]))}


def benchmark(root):
    expected={'afiro':-4.6475314286e2,'sc50a':-6.4575077059e1,'sc50b':-7e1}
    for name,target in expected.items():
        A,b,c,kinds,_=read_mps(Path(root)/f'{name}.mps')
        start=time.perf_counter()
        try:
            result=solve_lp(A,b,c,kinds)
            duration=time.perf_counter()-start
            print({'problem':name,'rows':len(b),'columns':len(c),'objective':result['objective'],
                   'netlib_reference':target,'abs_error':abs(result['objective']-target),
                   'max_residual':result['max_residual'],'phase1_iterations':result['phase1_iterations'],
                   'phase2_iterations':result['phase2_iterations'],'seconds':duration})
        except Exception as exc: print({'problem':name,'error':repr(exc)})
