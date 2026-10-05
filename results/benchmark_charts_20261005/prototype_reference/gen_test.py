import random,itertools,subprocess,json,sys,numpy as np
from scipy.optimize import minimize
def mps(name,c,Q,A,lo,up,ub,isint,sense='MIN'):
    n=len(c);m=len(A)
    L=[f'NAME {name}','ROWS',' N OBJ']+[f' L R{i}' for i in range(m)]+['COLUMNS']
    for j in range(n):
        if isint[j]: L.append(f"    M{j}  'MARKER'  'INTORG'")
        ents=[('OBJ',c[j])]+[(f'R{i}',A[i][j]) for i in range(m) if A[i][j]]
        for k in range(0,len(ents),2): L.append(f'    X{j} '+' '.join(f'{a} {b!r}' for a,b in ents[k:k+2]))
        if isint[j]: L.append(f"    M{j}e  'MARKER'  'INTEND'")
    L+=['RHS']+[f'    RHS R{i} {ub[i]!r}' for i in range(m)]
    L+=['BOUNDS']+[f' UP BND X{j} {up[j]!r}' for j in range(n)]
    L+=['QUADOBJ']+[f'    X{j} X{k} {Q[j][k]!r}' for j in range(n) for k in range(j+1)]
    L+=['ENDATA']
    return '\n'.join(L)+'\n'
def ref(c,Q,A,ub,up,isint):
    n=len(c);ints=[j for j in range(n) if isint[j]];best=None
    for combo in itertools.product(*[range(int(up[j])+1) for j in ints]):
        fixed=dict(zip(ints,combo));free=[j for j in range(n) if j not in fixed]
        x0=np.zeros(n)
        for j,v in fixed.items():x0[j]=v
        def f(z):
            x=x0.copy();x[free]=z;return float(np.array(c)@x+0.5*x@np.array(Q)@x)
        cons=[{'type':'ineq','fun':(lambda z,i=i:(ub[i]-np.array(A[i])@np.where(np.isin(range(n),free),0,x0)-np.array(A[i])[free]@z))} for i in range(len(A))]
        if not free:
            x=x0
            if all(np.array(A[i])@x<=ub[i]+1e-9 for i in range(len(A))): v=f([]);best=v if best is None or v<best else best
            continue
        r=minimize(f,np.zeros(len(free)),constraints=cons,bounds=[(0,up[j]) for j in free],method='SLSQP',options={'ftol':1e-12,'maxiter':500})
        if r.success:
            x=x0.copy();x[free]=r.x
            if all(np.array(A[i])@x<=ub[i]+1e-6 for i in range(len(A))):
                best=r.fun if best is None or r.fun<best else best
    return best
bad=0;tot=0
for seed in range(int(sys.argv[1]),int(sys.argv[2])):
    rnd=random.Random(26000+119+seed);n=rnd.randint(3,5);m=rnd.randint(1,3)
    B=[[rnd.uniform(-1,1) for _ in range(n)] for _ in range(n)]
    Q=[[round(sum(B[k][i]*B[k][j] for k in range(n)),6) for j in range(n)] for i in range(n)]
    for i in range(n):Q[i][i]+=0.05
    c=[round(rnd.uniform(-6,2),4) for _ in range(n)];A=[[round(rnd.uniform(0,2),3) for _ in range(n)] for _ in range(m)]
    up=[rnd.randint(2,4) for _ in range(n)];ub=[round(rnd.uniform(2,8),3) for _ in range(m)]
    isint=[rnd.random()<0.6 for _ in range(n)]
    if not any(isint):isint[0]=True
    open('m.mps','w').write(mps('t',c,Q,A,None,up,ub,isint))
    p=subprocess.run(['/tmp/r2/taral','m.mps','--time-limit','30','--json','o.json'],capture_output=True,text=True)
    j=json.load(open('o.json'));r=ref(c,Q,A,ub,up,isint);tot+=1
    got=j.get('objective') if j['status']=='optimal' else None
    ok=(r is None and j['status']=='infeasible') or (r is not None and got is not None and abs(got-r)<=1e-5*max(1,abs(r)))
    if not ok and got is not None and r is not None and got<r:
        x=np.zeros(n)
        for l in open('o.sol'):
            a=l.split()
            if len(a)==2 and a[0].startswith('X'):x[int(a[0][1:])]=float(a[1])
        fe=all(np.array(A[i])@x<=ub[i]+1e-6 for i in range(m)) and (x>=-1e-6).all() and (x<=np.array(up)+1e-6).all() and all(abs(x[k]-round(x[k]))<1e-6 for k in range(n) if isint[k])
        ov=float(np.array(c)@x+0.5*x@np.array(Q)@x)
        if fe and abs(ov-got)<1e-6: ok=True; better=globals().get('better',0)+1; globals()['better']=better
    if not ok: bad+=1;print('MISMATCH',seed,j['status'],got,r,p.stdout[:150])
print('tested',tot,'mismatch',bad,'ours_better_verified_feasible',globals().get('better',0))
