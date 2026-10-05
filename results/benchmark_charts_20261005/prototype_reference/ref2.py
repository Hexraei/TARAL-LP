import itertools,subprocess,json,sys,random,numpy as np
sys.argv=['x','0','1']
exec(open('gen_test.py').read().split('bad=0;tot=0')[0])
def fixed_mps(c,Q,A,ub,up,isint,fix):
    n=len(c);lo=[fix.get(j,0) for j in range(n)];hi=[fix.get(j,up[j]) for j in range(n)]
    s=mps('f',c,Q,A,None,hi,ub,[False]*n)
    if fix: s=s.replace('BOUNDS\n','BOUNDS\n'+''.join(f' LO BND X{j} {lo[j]!r}\n' for j in fix))
    return s
def solve_sub(c,Q,A,ub,up,isint,fix,tag):
    open(f'{tag}.mps','w').write(fixed_mps(c,Q,A,ub,up,isint,fix))
    subprocess.run(['./taral',f'{tag}.mps','--time-limit','20','--sol',f'{tag}.sol','--json',f'{tag}.json'],capture_output=True)
    j=json.load(open(f'{tag}.json'))
    if j['status']!='optimal': return j['status'],None,None
    n=len(c);x=np.zeros(n)
    for l in open(f'{tag}.sol'):
        a=l.split()
        if len(a)==2 and a[0].startswith('X'):x[int(a[0][1:])]=float(a[1])
    return 'optimal',x,None
def kkt_ok(c,Q,A,ub,up,fix,x):
    # independent verify: feasibility + optimality via projected gradient condition on free vars using scipy NNLS-free check: compare against dual by LP over gradient
    from scipy.optimize import linprog
    n=len(c);Qm=np.array(Q);g=np.array(c)+Qm@x;free=[j for j in range(n) if j not in fix]
    # convex QP optimal iff x minimises linearization g'y over the feasible set (first-order optimality)
    lo=[fix.get(j,0) for j in range(n)];hi=[fix.get(j,up[j]) for j in range(n)]
    r=linprog(g,A_ub=np.array(A),b_ub=np.array(ub),bounds=list(zip(lo,hi)),method='highs')
    if r.status!=0:return False
    return float(g@x-r.fun)<=1e-6*(1+abs(float(g@x)))
def exact_ref(c,Q,A,ub,up,isint,tag='r'):
    n=len(c);ints=[j for j in range(n) if isint[j]];best=None;allok=True
    for combo in itertools.product(*[range(int(up[j])+1) for j in ints]):
        fix=dict(zip(ints,combo));st,x,_=solve_sub(c,Q,A,ub,up,isint,fix,tag)
        if st!='optimal':continue
        if not kkt_ok(c,Q,A,ub,up,fix,x):allok=False;continue
        v=float(np.array(c)@x+0.5*x@np.array(Q)@x)
        if best is None or v<best:best=v
    return best,allok
if __name__=='__main__':
    a,b=int(sys.argv[1]) if False else 0,0
