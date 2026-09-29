from milp.branch_and_bound import branch_and_bound
from scipy.optimize import milp,Bounds,LinearConstraint
from itertools import product
import numpy as np,json,time
rng=np.random.default_rng(26000+119);counts={'pass':0,'fail':0,'infeasible':0};cases=[]
for seed in range(30):
 n=4;A=rng.integers(0,6,size=(3,n)).astype(float);ub=rng.integers(1,4,size=n);c=rng.integers(-8,2,size=n).astype(float)
 b=np.array([rng.integers(2,max(3,int(row@ub))) for row in A],float)
 AA=np.vstack([A,np.eye(n)]);bb=np.r_[b,ub];k=['L']*len(bb)
 brute=[]
 for t in product(*(range(u+1) for u in ub)):
  x=np.array(t,float)
  if np.all(AA@x<=bb+1e-9):brute.append((float(c@x),x))
 expected=min(brute,key=lambda t:t[0])[0]
 try:
  o=branch_and_bound(AA,bb,c,k,range(n),node_limit=10000)
  h=milp(c,integrality=np.ones(n),bounds=Bounds(np.zeros(n),ub),constraints=LinearConstraint(A,-np.inf*np.ones(3),b))
  passed=o['status']=='optimal' and h.success and abs(o['objective']-expected)<1e-7 and abs(h.fun-expected)<1e-7 and all(abs(v-round(v))<1e-7 for v in o['x'])
  if not passed:cases.append({'seed':seed,'status':o['status'],'obj':o['objective'],'brute':expected,'highs':h.fun})
  counts['pass' if passed else 'fail']+=1
 except Exception as ex:cases.append({'seed':seed,'error':repr(ex),'brute':expected});counts['fail']+=1
print(json.dumps({'tests':30,'counts':counts,'failures':cases}))
