#!/usr/bin/env python3
# R2 convex-MIQP harness: seeded generator + exact-enumeration reference + proto comparison.
# Corpus: random separable convex MIQPs (diagonal Q >= 0), feasible by construction.
# Seeds: 26119+i (program seed family 26000+119). Reference: exhaustive integer assignment
# enumeration; each assignment's continuous restriction solved by HiGHS QP (PSD diagonal Hessian).
# Resume-safe: skips case indices already in the ledger. Chunked: --start --count.
import argparse, csv, itertools, json, os, random, subprocess, sys, time
PROTO="/home/sandbox/taral/proto/build/taral"
BASE="/home/sandbox/taral/r2"
ap=argparse.ArgumentParser(); ap.add_argument("--start",type=int,required=True); ap.add_argument("--count",type=int,required=True)
a=ap.parse_args()
LED=os.path.join(BASE,"r2_ledger.csv")
FIELDS=["idx","seed","n_int","n_cont","n_rows","proto_status","proto_obj","proto_nodes","ref_obj","ref_assignments","rel_diff","verdict"]
done=set()
if os.path.exists(LED):
    for r in csv.DictReader(open(LED)):
        if r["idx"].isdigit(): done.add(int(r["idx"]))
led=open(LED,"a",newline=""); w=csv.DictWriter(led,fieldnames=FIELDS)
if not done: w.writeheader(); led.flush()

def gen(idx):
    rng=random.Random(26119+idx)
    ni=rng.randint(3,7); nc=rng.randint(2,6); m=rng.randint(3,8)
    int_hi=[rng.choice([1,1,1,3]) for _ in range(ni)]   # mostly binary, some 0..3
    cont_hi=[rng.choice([2.0,5.0,10.0]) for _ in range(nc)]
    x0=[rng.randint(0,h) for h in int_hi]+[rng.uniform(0,h) for h in cont_hi]
    rows=[]
    for _ in range(m):
        arow=[rng.uniform(-2,2) for _ in range(ni+nc)]
        act=sum(ai*xi for ai,xi in zip(arow,x0))
        rows.append((arow, act-rng.uniform(0.5,5), act+rng.uniform(0.5,5)))
    q=[rng.uniform(0,2) if rng.random()<0.8 else 0.0 for _ in range(ni+nc)]
    c=[rng.uniform(-3,3) for _ in range(ni+nc)]
    return ni,nc,int_hi,cont_hi,rows,q,c

def write_mps(path,ni,nc,int_hi,cont_hi,rows,q,c):
    L=["NAME          R2CASE","ROWS"," N  obj"]
    for i in range(len(rows)): L.append(f" G  r{i}")
    L+=["COLUMNS","    MARKER                 'MARKER'                 'INTORG'"]
    for j in range(ni):
        L.append(f"    x{j} obj {c[j]:.10g}")
        for i,(ar,lo,up) in enumerate(rows): L.append(f"    x{j} r{i} {ar[j]:.10g}")
    L+=["    MARKER                 'MARKER'                 'INTEND'"]
    for j in range(ni,ni+nc):
        L.append(f"    x{j} obj {c[j]:.10g}")
        for i,(ar,lo,up) in enumerate(rows): L.append(f"    x{j} r{i} {ar[j]:.10g}")
    L.append("RHS")
    for i,(ar,lo,up) in enumerate(rows): L.append(f"    RHS r{i} {lo:.10g}")
    L.append("RANGES")
    for i,(ar,lo,up) in enumerate(rows): L.append(f"    RNG r{i} {up-lo:.10g}")
    L.append("QUADOBJ")
    for j in range(ni+nc):
        if q[j]>0: L.append(f"    x{j} x{j} {2*q[j]:.10g}")
    L.append("BOUNDS")
    for j in range(ni): L.append(f" UP BND x{j} {int_hi[j]}")
    for j in range(ni,ni+nc): L.append(f" UP BND x{j} {cont_hi[j-ni]:.10g}")
    L.append("ENDATA")
    # G row with RANGES: G means a'x >= lo; range makes it <= lo+(up-lo). Standard MPS semantics.
    open(path,"w").write("\n".join(L)+"\n")

def enum_ref(ni,nc,int_hi,cont_hi,rows,q,c):
    import numpy as np
    from scipy.optimize import minimize, LinearConstraint, Bounds
    best=None; nassign=0
    doms=[range(h+1) for h in int_hi]
    A=np.array([ar for (ar,lo,up) in rows])
    for assign in itertools.product(*doms):
        nassign+=1
        f=A[:, :ni] @ np.array(assign,float)
        rl=np.array([lo for (ar,lo,up) in rows])-f
        ru=np.array([up for (ar,lo,up) in rows])-f
        Ac=A[:, ni:]
        def fun(x):
            return 0.5*np.sum(2*np.array(q[ni:])*x*x)+np.dot(c[ni:],x)
        def jac(x):
            return 2*np.array(q[ni:])*x+np.array(c[ni:])
        lc=LinearConstraint(Ac, rl, ru)
        r=minimize(fun, np.zeros(nc), jac=jac, bounds=Bounds(np.zeros(nc), np.array(cont_hi)),
                   constraints=[lc], method="SLSQP", options={"maxiter":200,"ftol":1e-12})
        if not r.success: continue
        obj=float(r.fun)+sum(q[j]*assign[j]**2 for j in range(ni))+sum(c[j]*assign[j] for j in range(ni))
        if best is None or obj<best: best=obj
    return best, nassign

for idx in range(a.start, a.start+a.count):
    if idx in done: continue
    ni,nc,int_hi,cont_hi,rows,q,c=gen(idx)
    mp=os.path.join(BASE,"cases",f"case{idx}.mps")
    write_mps(mp,ni,nc,int_hi,cont_hi,rows,q,c)
    jp=mp+".json"
    try:
        subprocess.run([PROTO,mp,"--time-limit","60","--json",jp],timeout=70,capture_output=True)
        j=json.load(open(jp))
    except Exception as e:
        j={"status":"runner_error","message":str(e)}
    robj,nas=enum_ref(ni,nc,int_hi,cont_hi,rows,q,c)
    ps=j.get("status"); po=j.get("objective") if j.get("has_solution") else None
    if robj is None:
        v="ref_infeasible_skip" if ps=="infeasible" else "REF_FAIL"
        rd=None
    elif ps=="optimal" and po is not None:
        rd=abs(po-robj)/max(1.0,abs(robj))
        v="agree" if rd<=1e-6 else "DISAGREE"
    else:
        v=f"proto_{ps}"; rd=None
    w.writerow({"idx":idx,"seed":26119+idx,"n_int":ni,"n_cont":nc,"n_rows":len(rows),
        "proto_status":ps,"proto_obj":po,"proto_nodes":j.get("nodes"),"ref_obj":robj,
        "ref_assignments":nas,"rel_diff":rd,"verdict":v}); led.flush()
    print(f"[{time.strftime('%H:%M:%S')}] case{idx}: {ps} vs ref -> {v}",flush=True)
print("CHUNK_COMPLETE",flush=True)
