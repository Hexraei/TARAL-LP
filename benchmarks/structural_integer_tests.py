#!/usr/bin/env python3
"""CLI parity proofs + negative/corrupted-proof tests. No benchmark-name dispatch."""
import sys,json,tempfile,subprocess,copy,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import structural_integer_replay as replay
from milp_tests import Case,write_mps
import numpy as np
BIN=sys.argv[1];fixtures=Path(sys.argv[2]);work=Path(tempfile.mkdtemp());passed=0

def run(p,flag=True):
    js=work/(p.stem+'.json')
    subprocess.run([BIN,str(p),'--time-limit','5','--json',str(js)]+(['--integer-structure'] if flag else []),capture_output=True,text=True,check=False,timeout=15)
    return json.load(open(js))
for name,st,obj in [('enlight8','optimal',27),('enlight_hard','optimal',37),('enlight9','infeasible',None),('enlight11','infeasible',None)]:
    p=fixtures/(name+'.mps');j=run(p);assert j['status']==st
    if obj is not None:assert j['objective']==obj
    M=replay.load(str(p));replay.check(M,j);passed+=1
    for mut in ['rows','type','point','status']:
        bad=copy.deepcopy(j);c=bad['structural_certificate']
        if mut=='rows':
            if 'rows' in c:c['rows']=[]
            else:c['fixed_binary'][0]['rows']=[]
        elif mut=='type':c['type']='unknown'
        elif mut=='point':
            if 'x' in bad:bad['x'][0]+=1
            else:c['rows']=[0]
        else:bad['status']='optimal' if st=='infeasible' else 'infeasible'
        try:replay.check(M,bad)
        except (AssertionError,KeyError):passed+=1
        else:raise AssertionError('corrupted proof accepted '+name+' '+mut)
# Unsupported examples must give the same answer and no structural certificate.
examples=[
 ('fractional',[1],[[.5]],[.5],[.5],[0],[1],[True]),
 ('general_residue',[1,0],[[1,-2]],[2],[2],[0,0],[math.inf,math.inf],[True,True]),
 ('continuous',[1,0],[[1,-2]],[1],[1],[0,0],[1,math.inf],[True,False]),
 ('aux_violated',[0,0],[[1,-2]],[-3],[-3],[0,0],[1,0],[True,True]),
 ('free_binary',[1,1],[[1,1]],[1],[1],[0,0],[1,1],[True,True]),
 ('decimal_near_int',[1],[[1.0000000001]],[1],[1],[0],[1],[True]),
]
for name,c,A,rl,ru,lo,up,it in examples:
    p=work/(name+'.mps');write_mps(Case(name,c,A,rl,ru,lo,up,it),str(p))
    a=run(p);b=run(p,False);assert 'structural_certificate' not in a,name
    assert a['status']==b['status'] and a['objective']==b['objective'];passed+=1
# General all-integer even equation with odd RHS must prove infeasibility, not require binaries.
p=work/'general_odd.mps';write_mps(Case('general_odd',[1,0],[[2,4]],[1],[1],[0,0],[10,10],[True,True]),str(p))
j=run(p);assert j['status']=='infeasible';replay.check(replay.load(str(p)),j);passed+=1
# A normal CLI deadline applies; unsupported presolve combination cannot launder row indices.
r=subprocess.run([BIN,str(fixtures/'enlight8.mps'),'--integer-structure','--presolve'],capture_output=True);assert r.returncode==2;passed+=1
print('PASS',passed,'structural cases and corruption guards')
