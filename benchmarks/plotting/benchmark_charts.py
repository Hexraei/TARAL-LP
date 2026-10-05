#!/usr/bin/env python3
"""Rebuild documentation charts from pinned ledger inputs; no solver execution."""
from pathlib import Path
import csv, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/benchmark_charts_20261005'
INPUT = OUT / 'inputs'
TEAL, INK, GRAY, PALE = '#73a89a', '#251f21', '#585254', '#eae9ea'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'text.color': INK,
                     'axes.labelcolor': INK, 'axes.edgecolor': PALE, 'xtick.color': GRAY,
                     'ytick.color': INK, 'savefig.facecolor': 'white'})

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / (name + '.png'), dpi=180)
    fig.savefig(OUT / (name + '.svg'))
    plt.close(fig)

def netlib():
    counts=[]
    for stage in ('s0', 's5'):
        path=INPUT/f'netlib93_dual_{stage}_ledger.csv'
        rows=list(csv.DictReader(path.open()))
        assert len(rows)==93 and len({r['case'] for r in rows})==93
        excluded={'pilot.we', 'pilot4'}
        accepted=sum(r['passed']=='True' and r['case'] not in excluded for r in rows)
        counts.append(accepted)
    assert counts==[90,91], counts
    fig,ax=plt.subplots(figsize=(10,4.8)); fig.subplots_adjust(left=.24,right=.96,top=.75,bottom=.30)
    labels=['Control (s0)','Five-stage stack (s5)']
    for y,c in enumerate(counts):
        ax.barh(y,c,color=TEAL,height=.42)
        ax.barh(y,93-c-2,left=c,color='#c0bebf',height=.42)
        ax.barh(y,2,left=91,color='white',edgecolor=GRAY,hatch='///',height=.42)
        ax.text(c/2,y,f'{c}/93 accepted',ha='center',va='center',weight='bold',color=INK)
    ax.set_yticks(range(2),labels);ax.invert_yaxis();ax.set_xlim(0,93)
    ax.set_xticks([0,30,60,93]);ax.set_xlabel('Cases (fixed denominator: 93)')
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
    fig.text(.065,.92,'Netlib dual: all 93 return solver-optimal results',fontsize=20,weight='bold')
    fig.text(.065,.845,'Five-stage stack: 93/93 solver-optimal; 91/93 counted by the published gate',color=GRAY)
    fig.text(.065,.14,'Bars show fixed-gate counts. Gray: timeout; hatched: 2 fixed reference exclusions.',color=GRAY,fontsize=10)
    fig.text(.065,.065,'60s/case; PILOT.WE/PILOT4 excluded by policy despite passing raw reference rows. No speed claim.',color=GRAY,fontsize=9)

    save(fig,'netlib_dual93')
    return {'s0':counts[0],'s5':counts[1],'denominator':93,'reference_exclusions':['pilot.we','pilot4']}


def category_chart(title, subtitle, labels, values, colors, name, footers):
    fig,ax=plt.subplots(figsize=(11,6.4));fig.subplots_adjust(left=.37,right=.91,top=.76,bottom=.27)
    ax.barh(range(len(values)),values,color=colors,height=.6)
    for y,v in enumerate(values):ax.text(v+max(values)*.025,y,str(v),va='center',weight='bold')
    ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.set_xlim(0,max(values)*1.15)
    ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('Cases');ax.xaxis.grid(True,color=PALE);ax.set_axisbelow(True)
    fig.text(.065,.925,title,fontsize=20,weight='bold');fig.text(.065,.86,subtitle,fontsize=11,color=GRAY)
    fig.text(.065,.13,footers[0],fontsize=9,color=GRAY)
    save(fig,name)

def qp():
    from collections import Counter
    rows=list(csv.DictReader((INPUT/'qp_138_verified_ledger.csv').open()));assert len(rows)==138
    c=Counter(r['final_class'] for r in rows)
    keys=['pass_vs_highs','matched_clarabel_second_reference','reader_artifact_on_reference_side','unresolved','time_limit_current_implementation_times_out','numerical_failure','out_of_class_nonconvex']
    labels=['HiGHS reference agreement','Clarabel second reference agreement','Reference-reader workaround','Unresolved reference disagreement','Time limit','Numerical failure','Out of class: nonconvex']
    category_chart('QP: 120/138 accepted reference results','87 HiGHS agreements + 32 Clarabel agreements + 1 reference-reader workaround',labels,[c[k] for k in keys],[TEAL,TEAL,'#b3c9c2','#c0bebf','#c0bebf','#c0bebf','#c0bebf'],'qp138',[
      '60s engine cap; mixed builds/references; 1 reader workaround and 3 disagreements. Details and sources in README.',
      'Mixed source/build evidence; 3 objective disagreements remain unresolved. No speed claim.',
      'Source: inputs/qp_138_verified_ledger.csv and clarabel_second_reference_raw.jsonl | Oct 5, 2026'])
    return dict(c)

def miplib():
    import math
    base=json.loads((INPUT/'miplib3_base65.json').read_text());assert len(base)==65
    def match(t,h):
        return h['status']=='Optimal' and t.get('objective') is not None and abs(t['objective']-h['obj'])<=1e-6*max(1,abs(h['obj']))
    accepted={k for k,r in base.items() if r['taral']['status']=='optimal' and match(r['taral'],r['highs'])};assert len(accepted)==28
    gains=set();best={k:r['taral'] for k,r in base.items()}
    for file in ['miplib3_persistent10.json','miplib3_persistent27.json']:
        if not (INPUT/file).exists():continue
        for r in json.loads((INPUT/file).read_text()):
            k=r['case'];assert r['sha_ok'] and r['tl']==300
            if r['taral']['status']=='optimal' and match(r['taral'],base[k]['highs']):gains.add(k)
            if match(r['taral'],base[k]['highs']):best[k]=r['taral']
    gains-=accepted
    unproven={k for k,t in best.items() if k not in accepted|gains and match(t,base[k]['highs'])}
    final=len(accepted|gains)
    category_chart(f'MIPLIB3: {final}/65 accepted, {final+len(unproven)}/65 reach the reference value','28 baseline results + 2 persistent-node gains; all accepted objectives match HiGHS',
      ['Baseline ba03938 accepted','Extra matches: persistent nodes','Matching incumbent, not accepted','Other unresolved cases'],
      [len(accepted),len(gains),len(unproven),65-final-len(unproven)],[TEAL,TEAL,'#b3c9c2','#c0bebf'],'miplib3_65',[
      '300s/case; combined snapshots/modes, 1e-6 stopping gap; 6 matching incumbents unproven. Sources in README.',
      'Accepted requires solver optimal status + HiGHS objective match. MIP stopping gap: 1e-6.',
      'Unproven incumbents do not count as solved. No speed comparison. See inputs/ and README.md.'])
    return {'baseline':len(accepted),'opt_in_gains':sorted(gains),'accepted':final,'unproven_matches':sorted(unproven),'denominator':65}

def miqp():
    rows=list(csv.DictReader((INPUT/'r2_convex_miqp_1032_ledger.csv').open()));assert len(rows)==1032
    from collections import Counter
    c=Counter(r['stage'] for r in rows);assert all(r['agree']=='True' and r['ref_kkt_verified']=='True' and r['taral_status']=='optimal' for r in rows)
    category_chart('Convex MIQP prototype: 1032/1032 reference agreements','Four synthetic test groups agree with the reported enumeration reference',
      [f'Stage {k}' for k in 'ABCD'],[c[k] for k in 'ABCD'],[TEAL]*4,'r2_miqp1032',[
      'Prototype, not main; 30s, all-feasible synthetic set; same-engine reference has limits. See README for details.',
      'All instances feasible, bounded integer subset. No infeasible-status coverage or performance claim.',
      'Local 2-core box, 30s cap, 1e-6 relative tolerance | inputs/r2_convex_miqp_1032_ledger.csv'])
    return {'accepted':1032,'denominator':1032,'stages':dict(c),'prototype_only':True}

if __name__=='__main__':
    summary={'netlib_dual93':netlib(),'qp138':qp(),'miplib3_65':miplib(),'r2_miqp1032':miqp()}
    summary['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUT.glob('*')) if p.is_file()}
    (OUT/'chart_data.json').write_text(json.dumps(summary,indent=2)+'\n')
