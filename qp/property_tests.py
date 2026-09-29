from qp_projected import solve_box_qp
from scipy.optimize import minimize,Bounds
import numpy as np,json
rng=np.random.default_rng(26000+119);failed=[];maxerr=0
for seed in range(30):
 n=4;M=rng.standard_normal((n,n));Q=M.T@M+np.eye(n);c=rng.uniform(-8,2,n);u=rng.uniform(.5,4,n)
 try:
  ours=solve_box_qp(Q,c,u,tol=1e-9,max_iter=100000)
  h=minimize(lambda x:.5*x@Q@x+c@x,np.zeros(n),jac=lambda x:Q@x+c,bounds=Bounds(np.zeros(n),u),method='SLSQP',options={'ftol':1e-12,'maxiter':1000})
  err=abs(ours['objective']-h.fun);maxerr=max(maxerr,err)
  if not h.success or err>1e-7 or ours['projected_gradient_inf']>1e-8:failed.append({'seed':seed,'objective_error':err,'reference_success':h.success,'grad':ours['projected_gradient_inf']})
 except Exception as e:failed.append({'seed':seed,'error':repr(e)})
print(json.dumps({'tests':30,'passes':30-len(failed),'failures':failed,'max_objective_error':maxerr}))
