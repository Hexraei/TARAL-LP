"""Convex box-QP, first-order projected-gradient prototype; no general linear constraints.
Min 0.5*x.T@Q@x + c@x with 0<=x<=upper and Q symmetric positive definite.
Step = 1/L exactly, L=max eigenvalue(Q), supports bound constraints only.
"""
import numpy as np

def solve_box_qp(Q,c,upper,tol=1e-10,max_iter=100000):
 Q=np.asarray(Q,float);c=np.asarray(c,float);upper=np.asarray(upper,float)
 n=len(c)
 if Q.shape!=(n,n) or upper.shape!=(n,) or np.any(upper<0):raise ValueError('dimension or bounds')
 if not np.allclose(Q,Q.T,atol=1e-12):raise ValueError('not symmetric')
 eig=np.linalg.eigvalsh(Q)
 if eig[0]<=0:raise ValueError('requires positive definite Q')
 x=np.zeros(n);step=1/eig[-1]
 for it in range(max_iter):
  y=np.clip(x-step*(Q@x+c),0,upper)
  projected_grad=np.max(np.abs(x-y))/step
  x=y
  if projected_grad<=tol:return {'x':x,'objective':float(.5*x@Q@x+c@x),'iterations':it+1,'projected_gradient_inf':float(projected_grad)}
 raise ArithmeticError(f'QP iteration limit: projected gradient={projected_grad}')
