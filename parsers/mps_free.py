from mps_fixed import read_mps_extended as _old
import numpy as np
from pathlib import Path

def read_mps_extended(path):
 # Replace FR variable declarations with two nonnegative columns in a
 # lossless fixed-field MPS transform. Both column coefficients appear under
 # signed cloned names; RHS and row semantics are unchanged.
 text=Path(path).read_text().splitlines()
 free={line[14:22].strip() for line in text if line[1:3]=='FR' and line[:1].isspace()}
 if not free:return _old(path)
 if len(free)!=len(set(free)):raise ValueError('duplicate FR names')
 output=[];section=''
 for line in text:
  if line and not line[0].isspace() and line.split()[0] in {'ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA'}:section=line.split()[0]
  if section=='COLUMNS' and line[4:12].strip() in free:
   original=line[4:12].strip()
   for suffix,sgn in [('P',1),('M',-1)]:
    copy=line[:4]+(original[:7]+suffix).ljust(8)+line[12:]
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
 return A,b,c,k,names,offset,lower
