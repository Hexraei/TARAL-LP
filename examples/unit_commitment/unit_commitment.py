#!/usr/bin/env python3
# Unit-commitment style dispatch MILP (textbook scale): G generators, T periods.
# Minimize running + startup cost s.t. demand met each period, capacity, min-up simplification.
# Writes MPS, solves with taral and HiGHS, compares. Deterministic (no randomness).
import subprocess, json, sys, os
G, T = 4, 12
cap   = [80, 120, 50, 200]      # MW max per generator
cvar  = [22.0, 18.0, 30.0, 12.0]# cost per MWh
cfix  = [40.0, 60.0, 25.0, 90.0]# cost per period when on
cstart= [300.0, 500.0, 150.0, 800.0]
demand= [150, 170, 190, 210, 240, 260, 280, 250, 220, 200, 180, 160]
lines, ng = [], G
def v(g,t,k): return f"{k}{g}_{t}"
lines.append("NAME          UC_DISPATCH\nOBJSENSE\n MINIMIZE\nROWS\n N  COST")
for t in range(T): lines.append(f" G  DEM{t}")
for g in range(G):
    for t in range(T):
        lines.append(f" L  CAP{g}_{t}"); lines.append(f" G  ST{g}_{t}")
lines.append("COLUMNS")
cols={}
def ent(c,r,vv): cols.setdefault(c,[]).append((r,vv))
for g in range(G):
    for t in range(T):
        p,u,s=v(g,t,'p'),v(g,t,'u'),v(g,t,'s')
        ent(p,"COST",cvar[g]); ent(p,f"DEM{t}",1.0); ent(p,f"CAP{g}_{t}",1.0)
        ent(u,"COST",cfix[g]); ent(u,f"CAP{g}_{t}",-cap[g]); ent(u,f"ST{g}_{t}",-1.0)
        ent(s,"COST",cstart[g]); ent(s,f"ST{g}_{t}",1.0)
        if t>0: ent(v(g,t-1,'u'),f"ST{g}_{t}",1.0)
for c,es in cols.items():
    for r,vv in es: lines.append(f"    {c} {r} {vv}")
lines.append("RHS")
for t in range(T): lines.append(f"    RHS       DEM{t}      {demand[t]}")
for g in range(G):
    for t in range(T): lines.append(f"    RHS       CAP{g}_{t}   0.0"); lines.append(f"    RHS       ST{g}_{t}   0.0")
lines.append("BOUNDS")
for g in range(G):
    for t in range(T):
        lines.append(f" UP BND       {v(g,t,'p'):10s}{cap[g]}")
        lines.append(f" BV BND       {v(g,t,'u')}"); lines.append(f" BV BND       {v(g,t,'s')}")
lines.append("ENDATA")
open("unit_commitment.mps","w").write("\n".join(lines)+"\n")
# taral
subprocess.run([os.environ.get("TARAL_ENGINE", "taral"),"unit_commitment.mps","--time-limit","300",
                "--json","taral_uc.json","--sol","taral_uc.sol"],check=True,capture_output=True)
j=json.load(open("taral_uc.json"))
# HiGHS
import highspy
h=highspy.Highs(); h.setOptionValue("output_flag",False); h.readModel("unit_commitment.mps"); h.run()
ho=h.getObjectiveValue(); hs=str(h.getModelStatus()).split(".")[-1]
d=abs(j["objective"]-ho)/max(1.0,abs(ho))
out={"model":{"generators":G,"periods":T,"rows":T+2*G*T,"cols":3*G*T},
     "taral":{"status":j["status"],"objective":j["objective"],"nodes":j["nodes"],"wall_s":j["wall_s"]},
     "highs":{"status":hs,"objective":ho},"rel_diff":d,
     "agree_1e-6": bool(d<=1e-6 and j["status"]=="optimal")}
json.dump(out,open("unit_commitment_result.json","w"),indent=1)
print(json.dumps(out,indent=1))
