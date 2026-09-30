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
 c=np.array([columns[x].get(obj,0.) for x in names],float)
 lower=np.array([bounds.get(x,[0,np.inf])[0] for x in names]);upper=np.array([bounds.get(x,[0,np.inf])[1] for x in names])
 if not np.all(np.isfinite(lower)):raise NotImplementedError('nonfinite lower variable bounds')
 if np.any(upper<lower):raise ValueError('inconsistent bounds')
 # Build row semantics first, preserving the original range-row order.
 rowmap={};bs=[];kinds=[]
 for r in rnames:
  kind=rows[r];value=rhs.get(r,0.);inds=[len(bs)];bs.append(value);kinds.append(kind)
  if r in ranges:
   t=ranges[r];v=abs(t);inds.append(len(bs))
   if kind=='L':bs.append(value-v);kinds.append('G')
   elif kind=='G':bs.append(value+v);kinds.append('L')
   elif t>=0:kinds[-1]='G';bs.append(value+v);kinds.append('L')
   else:kinds[-1]='L';bs.append(value-v);kinds.append('G')
  rowmap[r]=inds
 finite=np.where(np.isfinite(upper))[0];mcore=len(bs)
 # One dense allocation only; populate stored MPS entries, not m*n dict lookups.
 A=np.zeros((mcore+len(finite),len(names)),float)
 for j,name in enumerate(names):
  for r,value in columns[name].items():
   if r in rowmap:
    for i in rowmap[r]:A[i,j]=value
 b=np.array(bs+[upper[j]-lower[j] for j in finite],float)
 offset=float(c@lower);b[:mcore]-=A[:mcore]@lower
 if len(finite):A[mcore+np.arange(len(finite)),finite]=1
 kinds.extend(['L']*len(finite))
 return A,b,c,kinds,names,offset,lower
