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
        basic_set=set(basis); candidate=min((j for j in sorted(allowed) if j not in basic_set and reduced[j]<-tol),key=lambda j:(reduced[j],j),default=None)
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

"""Conservative exact structural presolve for nonnegative LPs. No solver imported."""
import numpy as np
_solve = solve_lp

# Phase 2 sparse path: CSC simplex matrix and sparse LU for basis systems.
# The MPS reader and conservative structural presolve remain unchanged.
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


_HARRIS={'on':True}


def _sparse_revised(A,b,c,basis,allowed=None,max_iter=10000,tol=1e-9):
 m,n=A.shape;basis=list(basis);allowed=set(range(n)) if allowed is None else set(allowed)
 if m==0:raise ValueError('empty basis')
 stall=0
 for it in range(max_iter):
  lu=splu(A[:,basis].tocsc());xB=lu.solve(b)
  if np.min(xB)<-max(1e-7,1e-12*np.max(np.abs(b))):raise ArithmeticError('lost primal feasibility')
  pi=lu.solve(c[basis],trans='T');reduced=c-np.asarray(A.T@pi).ravel()
  basic_set=set(basis); candidate=min((j for j in sorted(allowed) if j not in basic_set and reduced[j]<-tol),key=lambda j:(reduced[j],j),default=None)
  if candidate is None:
   x=np.zeros(n);x[basis]=xB
   return x,float(c@x),basis,it
  d=lu.solve(A[:,candidate].toarray().ravel())
  pivot_cutoff=max(1e-8,1e-7*np.max(np.abs(d)))
  pos=np.flatnonzero(d>pivot_cutoff)
  if len(pos)==0:raise ArithmeticError('unbounded')
  # Harris two-pass ratio test: pass 1 finds the max step allowed with a small
  # feasibility tolerance; pass 2 picks the largest pivot among rows within it.
  xp=np.maximum(xB[pos],0.);dp=d[pos];ftol=1e-9
  if not _HARRIS['on']:
   leave=int(min(((xp[t]/dp[t],basis[pos[t]],int(pos[t])) for t in range(len(pos))))[2])
  elif stall>=30:
   # Degenerate stall: fall back to the plain minimum-ratio rule (ties by basis index).
   leave=int(min(((xp[t]/dp[t],basis[pos[t]],int(pos[t])) for t in range(len(pos))))[2])
  else:
   step=np.min((xp+ftol)/dp)
   cand=pos[(xp/dp)<=step]
   leave=int(cand[np.argmax(d[cand])])
  stall=stall+1 if xp[list(pos).index(leave)]/d[leave]<=1e-12 else 0
  basis[leave]=candidate
 raise ArithmeticError('iteration limit')


def _sparse_core_lp(A,b,c,kinds):
 """Minimize c*x with L/G/E rows and x>=0; sparse two-phase simplex."""
 A=sparse.csc_matrix(A,dtype=float); b=np.asarray(b,float).copy();c=np.asarray(c,float)
 kinds=list(kinds);m,n=A.shape
 if b.shape!=(m,) or c.shape!=(n,) or len(kinds)!=m:raise ValueError('dimensions')
 if m==0:raise ValueError('empty rows unsupported')
 signs=np.where(b<0,-1.,1.);A=sparse.diags(signs,format='csc')@A;b*=signs
 for i in np.flatnonzero(signs<0):kinds[i]={'L':'G','G':'L','E':'E'}[kinds[i]]
 for k in kinds:
  if k not in {'L','G','E'}:raise ValueError(k)
 slack_rows=[i for i,k in enumerate(kinds) if k!='E']
 slacks=sparse.csc_matrix((np.array([1 if kinds[i]=='L' else -1 for i in slack_rows],float),(slack_rows,np.arange(len(slack_rows)))),shape=(m,len(slack_rows)))
 core=sparse.hstack([A,slacks],format='csc');astart=core.shape[1]
 basis=[None]*m
 for j,i in enumerate(slack_rows):
  if kinds[i]=='L':basis[i]=n+j
 art_rows=[i for i in range(m) if basis[i] is None]
 arts=sparse.csc_matrix((np.ones(len(art_rows)),(art_rows,np.arange(len(art_rows)))),shape=(m,len(art_rows)))
 for j,i in enumerate(art_rows):basis[i]=astart+j
 full=sparse.hstack([core,arts],format='csc')
 phase1=np.r_[np.zeros(astart),np.ones(len(art_rows))]
 x,penalty,basis,it1=_sparse_revised(full,b,phase1,basis)
 if penalty>1e-7:raise ArithmeticError(f'infeasible; phase-I penalty={penalty}')
 keep=list(range(m))
 for i in range(m-1,-1,-1):
  if basis[i]<astart:continue
  lu=splu(full[keep,:][:,basis].tocsc())
  for j in range(astart):
   if j in basis:continue
   d=lu.solve(full[keep,j].toarray().ravel())
   if abs(d[i])>1e-8:
    basis[i]=j;break
  if basis[i]>=astart:
   keep.pop(i);basis.pop(i)
 phase2=np.r_[c,np.zeros(astart-n)]
 x,z,basis,it2=_sparse_revised(core[keep,:],b[keep],phase2,basis)
 residual_vec=np.asarray(A@x[:n]).ravel()-b
 residual=float(max(0.,*(max(0.,v) if k=='L' else max(0.,-v) if k=='G' else abs(v) for v,k in zip(residual_vec,kinds))))
 return {'x':x[:n],'objective':z,'phase1_iterations':it1,'phase2_iterations':it2,'kept_rows':len(keep),'max_residual':residual,'matrix_nnz':int(A.nnz),'sparse_basis_lu':True}

_validated_dense_core = _solve

def _sparse_solve_lp(A,b,c,kinds):
 # Sparse LU can lose feasibility for an ill-conditioned basis. If and only if
 # that happens, retry with the former independently validated dense LU core.
 # Never mark the sparse failure itself a pass; the dense retry must solve.
 try:
  return _sparse_core_lp(A,b,c,kinds)
 except ArithmeticError as exc:
  if str(exc) not in ('lost primal feasibility','iteration limit'):
   raise
  # Retry the sparse path once with the plain minimum-ratio rule (the v19 behavior).
  _HARRIS['on']=False
  try:
   r=_sparse_core_lp(A,b,c,kinds)
   return dict(r,harris_retry=True)
  except ArithmeticError as exc2:
   exc=exc2
   if str(exc) != 'lost primal feasibility':
    raise
  finally:
   _HARRIS['on']=True
  result=_validated_dense_core(A,b,c,kinds)
  return dict(result,sparse_fallback=True,sparse_fallback_reason=str(exc))

_solve = _sparse_solve_lp

def solve_lp(A,b,c,kinds, *, tol=1e-12):
 A=np.asarray(A,float); b=np.asarray(b,float); c=np.asarray(c,float); kinds=list(kinds)
 m,n=A.shape
 if b.shape!=(m,) or c.shape!=(n,) or len(kinds)!=m:raise ValueError('dimensions')
 origA=A;origb=b;origk=kinds;c_original=c
 # A matrix copy: never change the caller's arrays.
 A=A.copy();b=b.copy();c=c.copy()
 # Zero rows are feasibility tests, not pivots.
 active=np.any(A!=0,axis=1)
 for i in np.where(~active)[0]:
  ok=(0<=b[i]+tol if kinds[i]=='L' else 0>=b[i]-tol if kinds[i]=='G' else abs(b[i])<=tol)
  if not ok:raise ArithmeticError('infeasible; empty row')
 A=A[active];b=b[active];k=[kinds[i] for i in np.where(active)[0]]
 lower=np.zeros(n);upper=np.full(n,np.inf)
 # Detect single-variable constraints and fixed variables with exact structural zeros.
 singleton=[]
 for i in range(len(b)):
  nz=np.flatnonzero(A[i]);
  if len(nz)!=1:continue
  singleton.append(i);j=nz[0];v=b[i]/A[i,j];kind=k[i]
  if A[i,j]<0:kind={'L':'G','G':'L','E':'E'}[kind]
  if kind in ('G','E'):lower[j]=max(lower[j],v)
  if kind in ('L','E'):upper[j]=min(upper[j],v)
 if np.any(lower>upper+tol):raise ArithmeticError('infeasible; singleton bounds')
 fixed=np.where(np.isfinite(upper)&(abs(upper-lower)<=tol))[0]
 # Fixed values use the stricter original lower when values differ within tolerance.
 vals=lower[fixed]; xbase=np.zeros(n);xbase[fixed]=vals
 keepcol=np.ones(n,bool);keepcol[fixed]=False
 objshift=float(c[fixed]@vals)
 if len(fixed):b=b-A[:,fixed]@vals
 A=A[:,keepcol];c=c[keepcol]
 # Remove newly empty rows, checking feasibility after substitution.
 active=np.any(A!=0,axis=1)
 for i in np.where(~active)[0]:
  ok=(0<=b[i]+tol if k[i]=='L' else 0>=b[i]-tol if k[i]=='G' else abs(b[i])<=tol)
  if not ok:raise ArithmeticError('infeasible; empty row after fixed vars')
 A=A[active];b=b[active];k=[k[i] for i in np.where(active)[0]]
 # Zero columns: if cost negative, LP unbounded; else select x=0.
 live=np.any(A!=0,axis=0) if len(b) else np.zeros(len(c),bool)
 if np.any((~live)&(c<-tol)):raise ArithmeticError('unbounded; empty negative-cost column')
 rem=np.flatnonzero(keepcol)[live]
 if len(b)==0:
  x=xbase; z=objshift; phases=(0,0);kept=0
 else:
  try:r=_solve(A[:,live],b,c[live],k)
  except np.linalg.LinAlgError:
   # Conservative fallback to the original formulation when presolve makes
   # an artificial cleanup basis singular; no false success can be inferred.
   r0=_solve(origA,origb,np.asarray(c_original,float),origk)
   return dict(r0, presolve_fixed=0, presolve_empty_rows=0, presolve_empty_cols=0, presolve_fallback=True)
  x=xbase.copy();x[rem]=r['x'];z=r['objective']+objshift;phases=(r['phase1_iterations'],r['phase2_iterations']);kept=r['kept_rows']
 residual=float(max([0.]+[(max(0.,origA[i]@x-origb[i]) if origk[i]=='L' else max(0.,origb[i]-origA[i]@x) if origk[i]=='G' else abs(origA[i]@x-origb[i])) for i in range(m)]))
 return {'x':x,'objective':z,'max_residual':residual,'phase1_iterations':phases[0],'phase2_iterations':phases[1],'kept_rows':kept,'presolve_fixed':len(fixed),'presolve_empty_rows':m-int(np.count_nonzero(active)),'presolve_empty_cols':int(np.count_nonzero(~live))}
