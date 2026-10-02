"""False-status tests: random small LPs with known status (feasible/optimal, infeasible, unbounded).
Checks the solver never returns a result for an infeasible or unbounded LP, and that returned optima match HiGHS.
Usage: python false_status_tests.py --engine engine.py [--n 300] [--seed 26000+119]"""
import argparse, importlib.util, numpy as np
from scipy.optimize import linprog
ap=argparse.ArgumentParser();ap.add_argument('--engine',required=True);ap.add_argument('--n',type=int,default=300);ap.add_argument('--seed',type=int,default=26000+119);a=ap.parse_args()
sp=importlib.util.spec_from_file_location('e',a.engine);E=importlib.util.module_from_spec(sp);sp.loader.exec_module(E)
rng=np.random.default_rng(a.seed);cnt={'optimal_match':0,'optimal_mismatch':0,'infeasible_caught':0,'infeasible_FALSE_OPTIMAL':0,'unbounded_caught':0,'unbounded_FALSE_OPTIMAL':0,'other_error_on_feasible':0,'skipped':0}
for t in range(a.n):
    m=int(rng.integers(2,9));n=int(rng.integers(m+1,14));A=rng.integers(-3,6,(m,n)).astype(float);kinds=list(rng.choice(['L','G','E'],m));x0=rng.integers(0,5,n).astype(float);b=A@x0
    for i,k in enumerate(kinds):
        if k=='L':b[i]+=rng.integers(0,3)
        if k=='G':b[i]-=rng.integers(0,3)
    c=rng.integers(-4,5,n).astype(float)
    mode=t%3
    if mode==1: # infeasible: add contradictory pair
        j=int(rng.integers(0,n));row=np.zeros(n);row[j]=1;A=np.vstack([A,row,row]);b=np.r_[b,1.,3.];kinds+=['G','L'] if False else ['L','G']
        # x_j<=1 and x_j>=3
    if mode==2: c=-np.abs(c)-1  # push toward unboundedness when feasible region unbounded
    s=np.array(kinds);Aub=np.vstack([A[s=='L'],-A[s=='G']]);bub=np.r_[b[s=='L'],-b[s=='G']]
    r=linprog(c,A_ub=Aub if len(bub) else None,b_ub=bub if len(bub) else None,A_eq=A[s=='E'] if (s=='E').any() else None,b_eq=b[s=='E'] if (s=='E').any() else None,bounds=(0,None),method='highs')
    try: o=E.solve_lp(A,b,c,kinds);ok=True
    except Exception as e: ok=False;err=repr(e)
    if r.status==0:
        if ok and abs(o['objective']-r.fun)<=1e-6*max(1,abs(r.fun)) and o['max_residual']<=1e-6: cnt['optimal_match']+=1
        elif ok: cnt['optimal_mismatch']+=1
        else: cnt['other_error_on_feasible']+=1
    elif r.status==2:
        cnt['infeasible_FALSE_OPTIMAL' if ok else 'infeasible_caught']+=1
    elif r.status==3:
        cnt['unbounded_FALSE_OPTIMAL' if ok else 'unbounded_caught']+=1
    else: cnt['skipped']+=1
print(cnt)
bad=cnt['optimal_mismatch']+cnt['infeasible_FALSE_OPTIMAL']+cnt['unbounded_FALSE_OPTIMAL']
print('FALSE_VERDICTS',bad)
