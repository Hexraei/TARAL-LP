"""Small pure-integer/binary MILP via independent B&B over local LU LP relaxations.
No cuts, warm starts, large-instance guarantees or GPU path. All x>=0.
"""
import numpy as np, math, heapq, time
from engine.solver_lu_relfeas import solve_lp

def branch_and_bound(A,b,c,kinds,integer_indices,node_limit=10000,integrality_tol=1e-7):
 A=np.asarray(A,float);b=np.asarray(b,float);c=np.asarray(c,float)
 if A.ndim!=2 or A.shape[1]!=len(c) or len(b)!=len(kinds) or len(b)!=len(A):raise ValueError('dimensions')
 idx=sorted(set(integer_indices));n=len(c)
 if not idx or any(j<0 or j>=n for j in idx):raise ValueError('integer indices')
 # Depth-first branch-and-bound, all node bounds included as ordinary LP rows.
 stack=[([],[],[])];best=math.inf;incumbent=None;explored=0;pruned_bound=0;infeasible=0
 while stack:
  if explored>=node_limit:return dict(status='node_limit',objective=best if incumbent is not None else None,x=incumbent,nodes=explored,open_nodes=len(stack))
  addedA,addedb,addedk=stack.pop();explored+=1
  aa=np.vstack([A]+addedA) if addedA else A
  bb=np.r_[b,addedb] if addedb else b
  kk=list(kinds)+addedk
  try:r=solve_lp(aa,bb,c,kk)
  except ArithmeticError as exc:
   if 'infeasible' in str(exc):infeasible+=1;continue
   raise # never reclassify numerical/unbounded failures as certified prune
  x=r['x'];z=r['objective']
  if r['max_residual']>1e-6:raise ArithmeticError('uncertified LP relaxation')
  if z>=best-1e-8:pruned_bound+=1;continue
  fractional=[j for j in idx if abs(x[j]-round(x[j]))>integrality_tol]
  if not fractional:
   if not all(abs(x[j]-round(x[j]))<=integrality_tol for j in idx):raise ArithmeticError('integer check')
   incumbent=x.copy();best=z;continue
  j=max(fractional,key=lambda j:(abs(x[j]-round(x[j])),-j))
  lo=math.floor(x[j]);hi=math.ceil(x[j]);left=np.zeros(n);left[j]=1.;right=np.zeros(n);right[j]=1.
  # Branch right first via LIFO to try smaller values before larger.
  stack.append((addedA+[right],addedb+[hi],addedk+['G']))
  stack.append((addedA+[left],addedb+[lo],addedk+['L']))
 return dict(status='optimal' if incumbent is not None else 'infeasible',objective=best if incumbent is not None else None,x=incumbent,nodes=explored,pruned_bound=pruned_bound,infeasible_nodes=infeasible,open_nodes=0)
