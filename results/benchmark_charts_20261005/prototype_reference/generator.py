import sys,random,json,subprocess,csv,numpy as np
NW=int(sys.argv[1]);WID=int(sys.argv[2])
import ref2
from ref2 import *
def gen(seed):
    rnd=random.Random(26000+119+seed)
    if seed<700:st='A';n=rnd.randint(3,5);m=rnd.randint(1,3)
    elif seed<1200:st='B';n=rnd.randint(7,10);m=rnd.randint(3,6)
    elif seed<2100:st='C';n=rnd.randint(8,10);m=rnd.randint(3,6)
    else:st='D';n=rnd.randint(9,11);m=rnd.randint(3,6)
    B=[[rnd.uniform(-1,1) for _ in range(n)] for _ in range(n)]
    Q=[[round(sum(B[k][i]*B[k][j] for k in range(n)),6) for j in range(n)] for i in range(n)]
    for i in range(n):Q[i][i]+=0.05
    c=[round(rnd.uniform(-6,2),4) for _ in range(n)];A=[[round(rnd.uniform(0,2),3) for _ in range(n)] for _ in range(m)]
    if st=='A':
        up=[rnd.randint(2,4) for _ in range(n)];ub=[round(rnd.uniform(2,8),3) for _ in range(m)]
        isint=[rnd.random()<0.6 for _ in range(n)]
    else:
        up=[rnd.randint(2,3) for _ in range(n)];ub=[round(rnd.uniform(2,8),3) for _ in range(m)]
        isint=[False]*n;k={'B':rnd.randint(3,4),'C':5,'D':6}[st]
        if st=='B':pass
        else:pass
    return st,n,m,c,Q,A,up,ub,isint,rnd
seeds=[s for s in list(range(0,700))+list(range(1000,1200))+list(range(2000,2100))+list(range(3000,3032)) if s%NW==WID]
out=open(f'led_{WID}.jsonl','w')
for seed in seeds:
    rnd=random.Random(26000+119+seed)
    if seed<700:
        st='A';n=rnd.randint(3,5);m=rnd.randint(1,3)
    elif seed<1200:
        st='B';n=rnd.randint(7,10);m=rnd.randint(3,6)
    elif seed<2100:
        st='C';n=rnd.randint(8,10);m=rnd.randint(3,6)
    else:
        st='D';n=rnd.randint(9,11);m=rnd.randint(3,6)
    B=[[rnd.uniform(-1,1) for _ in range(n)] for _ in range(n)]
    Q=[[round(sum(B[k][i]*B[k][j] for k in range(n)),6) for j in range(n)] for i in range(n)]
    for i in range(n):Q[i][i]+=0.05
    c=[round(rnd.uniform(-6,2),4) for _ in range(n)];A=[[round(rnd.uniform(0,2),3) for _ in range(n)] for _ in range(m)]
    if st=='A':
        up=[rnd.randint(2,4) for _ in range(n)];ub=[round(rnd.uniform(2,8),3) for _ in range(m)]
        isint=[rnd.random()<0.6 for _ in range(n)]
        if not any(isint):isint[0]=True
    else:
        up=[rnd.randint(2,3) for _ in range(n)];ub=[round(rnd.uniform(2,8),3) for _ in range(m)]
        isint=[False]*n
        for k in rnd.sample(range(n),{'B':rnd.randint(3,4),'C':rnd.randint(5,5),'D':rnd.randint(6,6)}[st]):isint[k]=True
        if not any(isint):isint[0]=True
    open(f'm{WID}.mps','w').write(mps('t',c,Q,A,None,up,ub,isint))
    subprocess.run(['./taral',f'm{WID}.mps','--time-limit','30','--json',f'o{WID}.json'],capture_output=True)
    j=json.load(open(f'o{WID}.json'));got=j.get('objective') if j['status']=='optimal' else None
    r,allok=ref2.exact_ref(c,Q,A,ub,up,isint,tag=f'r{WID}')
    ok=(r is None and j['status']=='infeasible') or (r is not None and got is not None and abs(got-r)<=1e-6*max(1,abs(r)))
    out.write(json.dumps(dict(seed=seed,stage=st,n_vars=n,n_rows=m,n_int=sum(isint),n_cont=n-sum(isint),taral_status=j['status'],taral_obj=got,ref_obj=r,ref_kkt_verified=allok,agree=bool(ok and allok),nodes=j.get('nodes')))+'\n');out.flush()
