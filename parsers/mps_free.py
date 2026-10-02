from mps_fixed import read_mps_extended as _old
import numpy as np
from pathlib import Path

def read_mps_extended(path):
 # Replace FR variable declarations with two nonnegative columns in a
 # lossless fixed-field MPS transform. Both column coefficients appear under
 # signed clone columns; RHS and row semantics are unchanged. Clone names are
 # generated unique (never derived from a name prefix), and the returned names
 # are '<original>|FR+' / '<original>|FR-' so solutions map back exactly.
 text=Path(path).read_text().splitlines()
 free={line[14:22].strip() for line in text if line[1:3]=='FR' and line[:1].isspace()}
 if not free:return _old(path)
 existing=set();section=''
 for line in text:
  if line and not line[0].isspace():section=line.split()[0]
  elif section=='COLUMNS' and line.strip():existing.add(line[4:12].strip())
 clone={};k=0
 for original in sorted(free):
  pair=[]
  for suffix in 'PM':
   while True:
    k+=1;name=f'~{k:06d}'
    if name not in existing:break
   pair.append(name)
  clone[original]=tuple(pair)
 output=[];section=''
 for line in text:
  if line and not line[0].isspace() and line.split()[0] in {'ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA'}:section=line.split()[0]
  if section=='COLUMNS' and line[4:12].strip() in free:
   original=line[4:12].strip()
   for name,sgn in zip(clone[original],(1,-1)):
    copy=line[:4]+name.ljust(8)+line[12:]
    if sgn<0:
     chars=list(copy)
     for a,b in [(24,36),(49,61)]:
      val=copy[a:b].strip()
      if val:chars[a:b]=list(f'{-float(val.replace("D","E")):>12g}')
     copy=''.join(chars)
    output.append(copy)
  elif section=='BOUNDS' and line[1:3]=='FR' and line[14:22].strip() in free:
   # default nonnegative lower bound applies independently to cloned cols.
   continue
  else:output.append(line)
 temp=Path('/tmp')/(Path(path).stem+'-free.mps');temp.write_text('\n'.join(output)+'\n')
 A,b,c,k,names,offset,lower=_old(temp)
 back={name:f'{original}|FR{sign}' for original,pair in clone.items() for name,sign in zip(pair,'+-')}
 names=[back.get(n,n) for n in names]
 return A,b,c,k,names,offset,lower
