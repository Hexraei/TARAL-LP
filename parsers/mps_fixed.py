from pathlib import Path
import numpy as np

def read_mps_extended(path):
 section=None;rows={};columns={};rhs={};ranges={};bounds={};obj=None
 for line in Path(path).read_text().splitlines():
  if not line.strip() or line.startswith('*'):continue
  tok=line.split()
  if tok[0] in {'NAME','ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA'} and not line[0].isspace():
   section=tok[0];continue
  if section=='ROWS':
   kind,name=line[1:3].strip(),line[4:12].strip()
   if kind=='N':
    if obj is None:obj=name
   else:rows[name]=kind
  elif section=='COLUMNS':
   name=line[4:12].strip()
   if name not in columns:columns[name]={}
   pairs=[(line[14:22].strip(),line[24:36].strip()),(line[39:47].strip(),line[49:61].strip())]
   for key,v in pairs:
    if not key:continue
    val=float(v.replace('D','E'))
    columns[name][key]=columns[name].get(key,0.)+val
  elif section in ('RHS','RANGES'):
   table=rhs if section=='RHS' else ranges
   pairs=[(line[14:22].strip(),line[24:36].strip()),(line[39:47].strip(),line[49:61].strip())]
   for key,v in pairs:
    if key:table[key]=float(v.replace('D','E'))
  elif section=='BOUNDS':
   typ,name=line[1:3].strip(),line[14:22].strip();v=line[24:36].strip();val=float(v.replace('D','E')) if v else 0.
   old=bounds.get(name,[0.,np.inf])
   if typ=='LO':old[0]=val
   elif typ=='UP':old[1]=val
   elif typ=='FX':old=[val,val]
   elif typ=='FR':old=[-np.inf,np.inf]
   elif typ=='MI':old[0]=-np.inf
   elif typ=='PL':old[1]=np.inf
   elif typ=='BV':old=[0,1]
   else:raise NotImplementedError(f'bound {typ}')
   bounds[name]=old
 names=list(columns);rnames=list(rows)
 A=np.array([[columns[x].get(r,0.) for x in names] for r in rnames],float)
 b=np.array([rhs.get(r,0.) for r in rnames],float)
 c=np.array([columns[x].get(obj,0.) for x in names],float)
 kinds=list(rows.values());lower=np.array([bounds.get(x,[0,np.inf])[0] for x in names]);upper=np.array([bounds.get(x,[0,np.inf])[1] for x in names])
 if not np.all(np.isfinite(lower)):raise NotImplementedError('nonfinite lower variable bounds')
 if np.any(upper<lower):raise ValueError('inconsistent bounds')
 if ranges:
  rows2=[];b2=[];k2=[]
  for i,r in enumerate(rnames):
   rows2.append(A[i]);b2.append(b[i]);k2.append(kinds[i])
   if r not in ranges:continue
   t=ranges[r];v=abs(t)
   if kinds[i]=='L':rows2.append(A[i]);b2.append(b[i]-v);k2.append('G')
   elif kinds[i]=='G':rows2.append(A[i]);b2.append(b[i]+v);k2.append('L')
   else:
    # Replace equality by the appropriate first inequality; add other side.
    if t>=0:k2[-1]='G';rows2.append(A[i]);b2.append(b[i]+v);k2.append('L')
    else:k2[-1]='L';rows2.append(A[i]);b2.append(b[i]-v);k2.append('G')
  A=np.array(rows2);b=np.array(b2);kinds=k2
 # Substitute x = y+lower, y>=0. Add finite y<=upper-lower constraints.
 offset=float(c@lower);b=b-A@lower
 finite=np.where(np.isfinite(upper))[0]
 for j in finite:
  row=np.zeros(len(c));row[j]=1
  A=np.vstack((A,row));b=np.r_[b,upper[j]-lower[j]];kinds.append('L')
 return A,b,c,kinds,names,offset,lower
