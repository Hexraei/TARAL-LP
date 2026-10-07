#!/usr/bin/env python3
# Shift scheduling MILP (textbook scale): cover hourly staffing demand over a day with
# integer crews on 8-hour shifts starting at allowed hours; minimize total crew cost.
# Deterministic. Writes MPS (columns emitted contiguously), solves taral + HiGHS, compares.
import subprocess, json, os
H=24
demand=[4,3,2,2,2,3,4,6,8,9,9,8,8,7,6,6,7,8,9,8,7,6,5,4]  # required staff per hour
starts=[0,4,6,8,10,14,16]                                  # allowed shift start hours
L=8                                                        # shift length
cost=[60,60,60,70,70,80,60]                                # crew cost per shift
cols={}
for i,s in enumerate(starts):
    ent=[]
    for h in range(s,min(s+L,H)): ent.append((f"DEM{h}",1.0))
    ent.append(("COST",float(cost[i])))
    cols[f"x{s}"]=ent
lines=["NAME          SHIFT_SCHED","OBJSENSE"," MINIMIZE","ROWS"," N  COST"]
for h in range(H): lines.append(f" G  DEM{h}")
lines.append("COLUMNS")
for c,es in cols.items():
    for r,v in es: lines.append(f"    {c} {r} {v}")
lines.append("RHS")
for h in range(H): lines.append(f"    RHS DEM{h} {demand[h]}")
lines.append("BOUNDS")
for s in starts: lines.append(f" UI BND x{s} 20"); lines.append(f" LI BND x{s} 0")
lines.append("ENDATA")
open("shift_scheduling.mps","w").write("\n".join(lines)+"\n")
# integer markers needed for UI/LI to be integer: use MARKER section
txt=open("shift_scheduling.mps").read()
txt=txt.replace("COLUMNS\n","COLUMNS\n    MARKER                 'MARKER'                 'INTORG'\n",1)
txt=txt.replace("RHS\n","    MARKER                 'MARKER'                 'INTEND'\nRHS\n",1)
open("shift_scheduling.mps","w").write(txt)
subprocess.run([os.environ.get("TARAL_ENGINE", "taral"),"shift_scheduling.mps","--time-limit","300",
                "--json","taral_ss.json","--sol","taral_ss.sol"],check=True,capture_output=True)
j=json.load(open("taral_ss.json"))
import highspy
h=highspy.Highs(); h.setOptionValue("output_flag",False); h.readModel("shift_scheduling.mps"); h.run()
ho=h.getObjectiveValue(); hs=str(h.getModelStatus()).split(".")[-1]
d=abs(j["objective"]-ho)/max(1.0,abs(ho))
out={"model":{"hours":H,"shift_starts":len(starts),"rows":H+1,"cols":len(starts)},
     "taral":{"status":j["status"],"objective":j["objective"],"nodes":j["nodes"],"wall_s":j["wall_s"]},
     "highs":{"status":hs,"objective":ho},"rel_diff":d,
     "agree_1e-6": bool(d<=1e-6 and j["status"]=="optimal")}
json.dump(out,open("shift_scheduling_result.json","w"),indent=1)
print(json.dumps(out,indent=1))
