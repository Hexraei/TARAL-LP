"""Independent original-model checker for TARAL-LP.
Parses the original MPS file with its own parser (no code shared with the solver's reader),
maps the solver's returned vector back to the original variables, and checks:
 row activity vs row sense/range, variable bounds, objective recomputed from original costs
 (including the objective-row RHS constant), vs the solver's reported objective.
Usage: python orig_check.py --engine PATH/TO/engine.py --corpus DIR [--out FILE] CASE [CASE...]
The engine file must define solve_lp(A,b,c,kinds). A combined notebook cell also defines
read_mps_extended(path); separated repo engines use the sibling parsers/mps_free.py.
Repo usage: python benchmarks/orig_check.py --engine engine/revised_simplex.py --corpus DIR afiro
Requires highspy (pip install highspy) for the reference solve of the ORIGINAL MPS file.
Pass rule per case: HiGHS status optimal on the original MPS, |checker objective - HiGHS objective| <= max(1e-6,1e-7*|HiGHS objective|),
max relative row violation <= 1e-6, max bound violation <= 1e-6. 60 s solver limit per case. Pass rule/tolerances are not configurable here.
Note: the solver-reported objective can differ from the checker objective by the objective-row RHS constant (reported separately as solver_objective)."""
import sys, os, json, time, multiprocessing as mp
from pathlib import Path
import numpy as np
INF=float('inf')
def parse_mps(path):
    try: return _parse_mps(path,False)
    except ValueError: return _parse_mps(path,True)
def _parse_mps(path,fixed):
    sec=None; rows={}; order=[]; obj=None; cols={}; corder=[]; rhs={}; rng={}; bnd={}
    for raw in Path(path).read_text().splitlines():
        if not raw.strip() or raw[0]=='*': continue
        t=raw.split()
        if raw[0]!=' ':
            sec=t[0]; continue
        if fixed:
            P=lambda a,b: raw[a:b].strip()
            if sec=='ROWS': t=[P(1,3),P(4,12)]
            elif sec=='BOUNDS': t=[P(1,3),P(4,12),P(14,22)]+([P(24,36)] if P(24,36) else [])
            elif sec in('COLUMNS','RHS','RANGES'):
                t=[P(4,12),P(14,22),P(24,36)]+([P(39,47),P(49,61)] if P(39,47) else [])
        f=lambda s: float(s.replace('D','E').replace('d','e'))
        if sec=='ROWS':
            if t[0]=='N':
                if obj is None: obj=t[1]
            else: rows[t[1]]=t[0]; order.append(t[1])
        elif sec=='COLUMNS':
            if t[1]=='\'MARKER\'': continue
            j=t[0]
            if j not in cols: cols[j]={}; corder.append(j)
            for k in range(1,len(t)-1,2): cols[j][t[k]]=cols[j].get(t[k],0.)+f(t[k+1])
        elif sec in('RHS','RANGES'):
            d=rhs if sec=='RHS' else rng
            rest=t[1:] if len(t)%2==1 else t
            for k in range(0,len(rest)-1,2): d[rest[k]]=f(rest[k+1])
        elif sec=='BOUNDS':
            typ=t[0]; name=t[2] if len(t)>=3 and t[1] not in cols and t[2] in cols else t[1]
            val=f(t[-1]) if len(t)>=4 or (len(t)==3 and t[2] not in cols) else 0.
            lo,up=bnd.get(name,[0.,INF])
            if typ=='LO': lo=val
            elif typ=='UP':
                up=val
                if val<0 and lo==0.: lo=-INF
            elif typ=='FX': lo=up=val
            elif typ=='FR': lo,up=-INF,INF
            elif typ=='MI': lo=-INF
            elif typ=='PL': up=INF
            elif typ=='BV': lo,up=0.,1.
            bnd[name]=[lo,up]
    return dict(rows=rows,order=order,obj=obj,cols=cols,corder=corder,rhs=rhs,rng=rng,bnd=bnd)
def check(model,xmap):
    rows=model['rows']; act={r:0. for r in model['order']}; cobj=0.
    x={j:xmap.get(j,0.) for j in model['corder']}
    for j,col in model['cols'].items():
        v=x[j]
        for r,a in col.items():
            if r==model['obj']: cobj+=a*v
            elif r in act: act[r]+=a*v
    const=-model['rhs'].get(model['obj'],0.)
    rowviol=0.; worst=None
    for r in model['order']:
        k=rows[r]; b=model['rhs'].get(r,0.); R=model['rng'].get(r)
        lo,up=-INF,INF
        if k=='L': up=b; lo=b-abs(R) if R is not None else -INF
        elif k=='G': lo=b; up=b+abs(R) if R is not None else INF
        else:
            if R is None: lo=up=b
            elif R>=0: lo,up=b,b+R
            else: lo,up=b+R,b
        v=max(lo-act[r],act[r]-up,0.)/(1.+abs(b))
        if v>rowviol: rowviol=v; worst=r
    bviol=0.
    for j in model['corder']:
        lo,up=model['bnd'].get(j,[0.,INF]); v=x[j]
        bviol=max(bviol,(lo-v) if v<lo else 0.,(v-up) if v>up else 0.)
    return dict(objective=cobj+const,row_violation_rel=rowviol,worst_row=worst,bound_violation=bviol)
def solve_and_check(path,E,q):
    try:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            from contextlib import nullcontext
            threadpool_limits=lambda limits: nullcontext()
        with threadpool_limits(limits=1):
            A,b,c,k,names,offset,lower=E.read_mps_extended(path)
            t=time.perf_counter(); r=E.solve_lp(A,b,c,k); wall=time.perf_counter()-t
        xv=np.asarray(r['x'],float)+np.asarray(lower,float)
        xmap={}
        for n,v in zip(names,xv): xmap[n]=xmap.get(n,0.)+v
        model=parse_mps(path)
        # free-variable clones: solver names are first 7 chars + P/M suffix
        for j in model['corder']:
            if j not in xmap and (j[:7]+'P') in xmap:
                xmap[j]=xmap.get(j[:7]+'P',0.)-xmap.get(j[:7]+'M',0.)
        res=check(model,xmap)
        res.update(solver_objective=float(r['objective']+offset),solver_residual=float(r['max_residual']),wall_s=wall)
        q.put(res)
    except BaseException as e:
        q.put({'error':repr(e)})
def highs_original(path):
    import highspy
    h=highspy.Highs(); h.setOptionValue('output_flag',False); h.setOptionValue('time_limit',60.0)
    h.readModel(str(path)); h.run()
    return str(h.getModelStatus()), h.getInfo().objective_function_value
def main(argv):
    import argparse, importlib.util
    ap=argparse.ArgumentParser(); ap.add_argument('--engine',required=True); ap.add_argument('--corpus',required=True)
    ap.add_argument('--out',default='orig_check_out.json'); ap.add_argument('cases',nargs='+'); a=ap.parse_args(argv)
    spec=importlib.util.spec_from_file_location('engine_under_test',a.engine); E=importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
    if not hasattr(E,'read_mps_extended'):
        # Repo engines keep the notebook parser in a separate module.
        parser_dir=Path(a.engine).resolve().parent.parent/'parsers'
        sys.path.insert(0,str(parser_dir))
        from mps_free import read_mps_extended
        E.read_mps_extended=read_mps_extended
    ctx=mp.get_context('fork'); out={}
    for n in a.cases:
        path=Path(a.corpus)/(n+'.mps'); q=ctx.Queue(1); p=ctx.Process(target=solve_and_check,args=(str(path),E,q)); p.start()
        try: r=q.get(timeout=60)
        except Exception: r={'status':'timeout'}
        p.join(2)
        if p.is_alive(): p.kill()
        if 'objective' in r:
            st,ho=highs_original(path); tol=max(1e-6,1e-7*abs(ho))
            r.update(highs_original_status=st,highs_original_objective=ho,
                     independent_pass=bool('ptimal' in st and abs(r['objective']-ho)<=tol and r['row_violation_rel']<=1e-6 and r['bound_violation']<=1e-6))
        else: r['independent_pass']=False
        out[n]=r; print(n,'PASS' if r['independent_pass'] else 'FAIL',json.dumps(r),flush=True)
    json.dump(out,open(a.out,'w'),indent=1)
    print('INDEPENDENT_PASS',sum(v['independent_pass'] for v in out.values()),'of',len(out))
if __name__=='__main__': main(sys.argv[1:])
