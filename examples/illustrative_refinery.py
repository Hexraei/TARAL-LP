"""Illustrative refinery planning LP, independently authored, not refinery operational data.
Units: thousand barrels per day (kbd), USD per barrel; quality scores are
synthetic indices. All coefficients and limits are assumptions, not refinery facts.
"""
import numpy as np
from scipy.optimize import linprog
import time, json
from engine.revised_simplex import solve_lp
# Crude L/H, each split into gasoline/diesel/residue by fixed yields.
crudes=['Light','Heavy'];products=['Gasoline','Diesel']
yields=np.array([[.42,.35,.18],[.27,.43,.25]]) # leftover 5% loss for each crude
crude_cost=np.array([72.,64.]);product_price=np.array([102.,95.,35.])
# Decision x_L,x_H,g_L,g_H,d_L,d_H; byproduct residue is deterministic.
# g_i,d_i <= crude yields; unused gasoline/diesel stream sold as low-value residue.
# More explicit independent blending streams with demands, quality, and line capacity.
names=['L','H','gL','gH','dL','dH'];n=len(names)
A=[];b=[];kind=[];labels=[]
def row(label,co,op,rhs):A.append(co);b.append(rhs);kind.append(op);labels.append(label)
row('CDU throughput <= 100 kbd',[1,1,0,0,0,0],'L',100)
row('Light supply <= 70 kbd',[1,0,0,0,0,0],'L',70)
row('Heavy supply <= 70 kbd',[0,1,0,0,0,0],'L',70)
row('Gasoline unit <= 40 kbd',[0,0,1,1,0,0],'L',40)
row('Diesel unit <= 42 kbd',[0,0,0,0,1,1],'L',42)
row('Gasoline demand >= 33 kbd',[0,0,1,1,0,0],'G',33)
row('Diesel demand >= 35 kbd',[0,0,0,0,1,1],'G',35)
for i,tag in enumerate(crudes):
 v=[0.]*n;v[i]=-yields[i,0];v[2+i]=1;row(f'{tag} gasoline yield',v,'L',0)
 v=[0.]*n;v[i]=-yields[i,1];v[4+i]=1;row(f'{tag} diesel yield',v,'L',0)
 # Product streams cannot in total exceed recoverable portions. Individually bounded above.
# Light blend RON 94, Heavy RON 88; minimum gasoline blend RON 91.
row('Gasoline RON >= 91',[0,0,94-91,88-91,0,0],'G',0)
# Light diesel sulfur 8 ppm, Heavy diesel sulfur 18 ppm; weighted average <= 14 ppm.
row('Diesel sulfur <= 14 ppm',[0,0,0,0,8-14,18-14],'L',0)
# Simplified mass balance: fixed residue 18/25% crude plus 5% loss and
# potentially unused gasoline/diesel yield. Unused clean streams are explicitly
# tracked as low-value byproduct sales in the linear objective below.
# No shipment/distribution/transport/energy costs are modeled.
# max margin: full recoverable yields (.95 per crude) sold as selected
# gasoline/diesel or at $35/bbl byproduct price; 5% per crude is lost.
# Revenue = 35*.95*(L+H) + (102-35)*(gL+gH) + (95-35)*(dL+dH).
c=np.array([72-35*.95,64-35*.95,-(102-35),-(102-35),-(95-35),-(95-35)],float)
A=np.array(A,float);b=np.array(b,float)
if __name__=='__main__':
 t=time.process_time();ours=solve_lp(A,b,c,kind);ot=time.process_time()-t
 le=np.array(kind)=='L';ge=np.array(kind)=='G'
 au=np.vstack([A[le],-A[ge]]);bu=np.r_[b[le],-b[ge]]
 t=time.process_time();hi=linprog(c,A_ub=au,b_ub=bu,bounds=(0,None),method='highs');ht=time.process_time()-t
 x=ours['x'];xh=hi.x
 def results(v):
  L,H,gL,gH,dL,dH=v;return {'decision_kbd':dict(zip(names,map(float,v))),
   'gross_margin_usd_per_day':float(-c@v*1000),
   'gasoline_kbd':float(gL+gH),'diesel_kbd':float(dL+dH),
   'fixed_residue_kbd':float(.18*L+.25*H),
   'total_byproduct_kbd':float(.95*(L+H)-(gL+gH+dL+dH)),
   'gasoline_ron':float((94*gL+88*gH)/(gL+gH)),
   'diesel_sulfur_ppm':float((8*dL+18*dH)/(dL+dH)),
   'max_constraint_violation_kbd':float(max([0]+[max(0,((a@v-rhs) if op=='L' else (rhs-a@v))) for a,rhs,op in zip(A,b,kind)]))}
 print(json.dumps({'case':'Illustrative refinery-flavored refinery LP, wholly synthetic',
  'variables':names,'objective_min_cost_per_kbd':c.tolist(),
  'constraints':[{'label':l,'coefficients':a.tolist(),'sense':k,'rhs':float(t)} for l,a,k,t in zip(labels,A,kind,b)],
  'crude_yields':yields.tolist(),'prices_usd_per_bbl':product_price.tolist(),
  'crude_cost_usd_per_bbl':crude_cost.tolist(),'ours':results(x),'highs':results(xh),
  'objective_difference_usd_per_day':float(abs(c@x-c@xh)*1000),
  'max_decision_difference_kbd':float(np.max(np.abs(x-xh))),
  'ours_cpu_s':ot,'highs_cpu_s':ht,'highs_status':int(hi.status),
  'iterations':{'phase1':ours['phase1_iterations'],'phase2':ours['phase2_iterations']}},indent=2))
