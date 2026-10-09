from fractions import Fraction as F
from pathlib import Path
import json,hashlib
D=Path(__file__).parent;p=D/'VALUES.mps';c=json.load(open(D/'exact_certificate.json'));assert hashlib.sha256(p.read_bytes()).hexdigest()==c['model_sha'];sec=None;names=[];row={};quad=[]
for line in p.read_text().splitlines():
 if not line.strip():continue
 ts=line.split()
 if not line[0].isspace():sec=ts[0];continue
 if sec=='COLUMNS':
  if ts[0]not in names:names.append(ts[0])
  for k in range(1,len(ts),2):
   if ts[k]=='r1':row[ts[0]]=F(ts[k+1])
 elif sec=='QUADOBJ':quad.append((ts[0],ts[1],F(ts[2])))
assert names==c['column_order'];d=dict(zip(names,c['direction_integer']));assert sum(row[n]*d[n]for n in names)==0
q=sum(v*d[a]*d[b]*(1 if a==b else 2)for a,b,v in quad);assert q==F(c['dQd_exact'])and q<0;t=F(c['step_exact'])
for sign in [-1,1]:
 x={n:F(5)+sign*t*d[n]for n in names};assert all(F(0)<=v<=F(10)for v in x.values());assert sum(row[n]*x[n]for n in names)==0
assert t*t*q==F(c['fplus_fminus_minus_2fmid_exact'])and t*t*q<0
print('PASS exact decimal-matrix negative-curvature witness and two feasible endpoints; no eigensolver/HiGHS dependency')
