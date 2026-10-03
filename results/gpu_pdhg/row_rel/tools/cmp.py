import csv,sys
old={r['case']:r for r in csv.DictReader(open(sys.argv[1]))}
new={r['case']:r for r in csv.DictReader(open(sys.argv[2]))}
f=lambda x:float(x) if x not in ('',None) else None
bad=0;ch=0
print(len(old),len(new))
for c in old:
    o,n=old[c],new.get(c)
    if n is None: print('MISSING',c);bad+=1;continue
    flags=[]
    if o['status']!=n['status']: flags.append(f"status {o['status']}->{n['status']}")
    eo,en=f(o['objective_rel_error_vs_highs']),f(n['objective_rel_error_vs_highs'])
    if eo is not None and en is not None and en>max(2*eo,1e-4 if '1e-4' in sys.argv[2] else 1e-6): flags.append(f"objerr {eo:.1e}->{en:.1e}")
    for k in ('chk_error',):
        if n.get(k): flags.append('CHKERR '+n[k])
    rv=f(n.get('chk_row_viol_rel')); 
    if n['iterations']!=o['iterations']: ch+=1
    if flags: bad+=1;print(c,flags)
print('rows with different iterations:',ch,' flagged:',bad)
print('max chk_row_viol_rel old/new:',max(f(r.get('chk_row_viol_rel')) or 0 for r in old.values()),max(f(r.get('chk_row_viol_rel')) or 0 for r in new.values()))
print('max chk_bound_viol new:',max(f(r.get('chk_bound_viol')) or 0 for r in new.values()))
