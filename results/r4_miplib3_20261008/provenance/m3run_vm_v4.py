#!/usr/bin/env python3
# MIPLIB3 65-case runner, Azure VM edition. Clean single-build ledger at c4800e9.
# Usage: python3 m3run_vm.py --shard 0 --shards 1   (shard i runs names[i::n])
# Protocol: taral --time-limit 300 --json --sol ; HiGHS reference 1 thread, 300 s. Resume-safe.
import json, subprocess, hashlib, os, csv, time, argparse
ap=argparse.ArgumentParser(); ap.add_argument("--shard",type=int,default=0); ap.add_argument("--shards",type=int,default=1)
ap.add_argument("--base",default=os.path.expanduser("~/r4run")); a=ap.parse_args()
TARAL=os.path.expanduser("~/r4run/TARAL-LP/build/taral"); BASE=a.base
LED=os.path.join(BASE,f"ledger_c4800e9_s{a.shard}.csv"); REF=os.path.join(BASE,f"ref_highs_s{a.shard}.json")
FIELDS=["instance","sha256","taral_status","taral_obj","best_bound","nodes","taral_wall_s",
        "highs_status","highs_obj","highs_wall_s","rel_diff","sol_max_viol","sol_feasible","verdict"]
os.makedirs(os.path.join(BASE,"out"),exist_ok=True)
done=set()
if os.path.exists(LED):
    for r in csv.DictReader(open(LED)): done.add(r["instance"])
led=open(LED,"a",newline=""); w=csv.DictWriter(led,fieldnames=FIELDS)
if not done: w.writeheader(); led.flush()
ref=json.load(open(REF)) if os.path.exists(REF) else {}
def highs_ref(name,path):
    if name in ref: return ref[name]
    import highspy
    h=highspy.Highs(); h.setOptionValue("output_flag",False)
    h.setOptionValue("threads",1); h.setOptionValue("time_limit",300.0); h.setOptionValue("mip_rel_gap",0.0)  # v4: exact-gap reference (misc06 false-WRONG fix)
    t=time.time(); h.readModel(path); h.run(); wall=time.time()-t
    st=str(h.getModelStatus()).split(".")[-1]
    obj=h.getObjectiveValue() if h.getSolution().value_valid else None
    ref[name]={"highs_status":st,"highs_obj":obj,"highs_wall_s":round(wall,3)}
    json.dump(ref,open(REF,"w")); return ref[name]
def sol_feasibility(mps_path, sol_path):
    # Independent point check: recompute original-row activities and bound violations from the
    # emitted .sol against the MPS as read by HiGHS (reader-independent of the engine's parser).
    try:
        import highspy
        h=highspy.Highs(); h.setOptionValue("output_flag",False); h.readModel(mps_path)
        lp=h.getLp()
        vals={}
        for ln in open(sol_path):
            p=ln.split()
            if len(p)==2: vals[p[0]]=float(p[1])
        names=h.getVariables()
        col_names=list(lp.col_names_) if lp.col_names_ else list(vals.keys())
        n=lp.num_col_
        x=[vals.get(col_names[j] if j<len(col_names) else "",0.0) for j in range(n)]
        import math
        act=[0.0]*lp.num_row_
        starts=list(lp.a_matrix_.start_); idx=list(lp.a_matrix_.index_); val=list(lp.a_matrix_.value_)
        for j in range(n):
            for k in range(starts[j],starts[j+1]): act[idx[k]]+=val[k]*x[j]
        worst=0.0
        for j in range(n):
            lo,up=lp.col_lower_[j],lp.col_upper_[j]
            if math.isfinite(lo): worst=max(worst,(lo-x[j])/(1+abs(lo)))
            if math.isfinite(up): worst=max(worst,(x[j]-up)/(1+abs(up)))
            integ=getattr(lp,"integrality_",None)
            if integ and integ[j]: worst=max(worst,abs(x[j]-round(x[j])))
        for i in range(lp.num_row_):
            lo,up=lp.row_lower_[i],lp.row_upper_[i]
            if math.isfinite(lo): worst=max(worst,(lo-act[i])/(1+abs(lo)))
            if math.isfinite(up): worst=max(worst,(act[i]-up)/(1+abs(up)))
        return worst, bool(worst<=1e-6)
    except Exception as e:
        return None, "check_error:"+str(e)[:60]

def classify(ts,tobj,hs,hobj,has_sol):
    if ts=="optimal" and hs=="kOptimal":
        d=abs(tobj-hobj)/max(1.0,abs(hobj))
        return ("proven_agree",d) if d<=1e-6 else ("WRONG_DISAGREE",d)
    if ts=="optimal": return ("proven_ref_unproven",None)
    if has_sol and hs=="kOptimal":
        d=abs(tobj-hobj)/max(1.0,abs(hobj))
        return ("value_reached_unproven",d) if d<=1e-6 else ("incumbent_differs_from_ref",d)
    return ("other",None)
names=[l.strip() for l in open(os.path.join(BASE,"list.txt")) if l.strip()]
mine=names[a.shard::a.shards]
for name in mine:
    if name in done: continue
    path=os.path.join(BASE,"mps",name+".mps")
    sha=hashlib.sha256(open(path,"rb").read()).hexdigest()
    jp=os.path.join(BASE,"out",name+".json"); sp=os.path.join(BASE,"out",name+".sol")
    try:
        subprocess.run([TARAL,path,"--time-limit","300","--json",jp,"--sol",sp],timeout=330,capture_output=True)
        j=json.load(open(jp))
    except Exception as e:
        j={"status":"runner_error","message":str(e)}
    r=highs_ref(name,path)
    sv,sf=(None,None)
    if j.get("has_solution") and os.path.exists(sp):
        sv,sf=sol_feasibility(path,sp)
    v,d=classify(j.get("status"),j.get("objective"),r["highs_status"],r["highs_obj"],j.get("has_solution",False))
    w.writerow({"instance":name,"sha256":sha,"taral_status":j.get("status"),"taral_obj":j.get("objective"),
        "best_bound":j.get("best_bound"),"nodes":j.get("nodes"),"taral_wall_s":j.get("wall_s"),
        "highs_status":r["highs_status"],"highs_obj":r["highs_obj"],"highs_wall_s":r["highs_wall_s"],
        "rel_diff":d,"sol_max_viol":sv,"sol_feasible":sf,"verdict":v}); led.flush()
    print(f"[{time.strftime('%H:%M:%S')}] s{a.shard} {name}: {j.get('status')} / ref {r['highs_status']} -> {v}",flush=True)
print("SHARD_COMPLETE",flush=True)
