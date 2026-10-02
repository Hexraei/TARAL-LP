#!/usr/bin/env python3
"""Round-3 PDHG harness for a fresh cloud GPU box (Modal, Colab; T4 or A10). Measured numbers only, near-optimal method.

Runs the rows of gpu/ROUND3_CASES.md (the single source of truth) with gpu/bench_pdhg.py and writes ledgers in the
schema of results/gpu_pdhg/synth_*.csv, one pair per family and tolerance (synth_<family>_tol<tol>.csv/json), plus
NOTES.md with the machine metadata, per-row outcome, failures, skips and deviations.

Builds only what PDHG needs: nvcc on gpu/pdhg.cu + src/mps.cpp. No CPU simplex, no HiGHS (references for known
cases are taken from the round-2 ledgers in results/gpu_pdhg/; new sizes have none).

Needs: nvcc, nvidia-smi, python3 with numpy, run from a checkout of this repository.
Usage: python3 gpu/run_round3.py [--work out/round3_work] [--out out/round3] [--budget-s 6600] [--only id,id] [--dry-run]
Resumable: a row whose rows/<id>.csv and .json exist is finished and skipped on re-run; failed and skipped rows are retried.
"""
import argparse, csv, json, os, platform, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
CASES_MD = os.path.join(HERE, 'ROUND3_CASES.md')
ROUND2 = os.path.join(ROOT, 'results', 'gpu_pdhg')
CHECK_RAM_GB = 24      # vector check of a 2M+ column model needs more host RAM than a small cloud box has
CHECK_MODEL_BYTES = 4.0e8


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def parse_cases(path):
    rows = []
    for line in open(path):
        c = [x.strip() for x in line.strip().strip('|').split('|')]
        if not line.startswith('|') or len(c) != 10 or c[0] == 'id' or set(c[0]) <= set('-'):
            continue
        rows.append(dict(id=c[0], family=c[1], size=c[2], tol=c[3], configs=c[4], cap=float(c[5]), repeats=int(c[6]),
                         check=c[7] == 'yes', est=float(c[8]), why=c[9]))
    if not rows:
        raise SystemExit('no case rows parsed from ' + path)
    return rows


def case_name(r):
    return f"{r['family']}_{r['size']}"


def machine():
    m = dict(python=platform.python_version(), platform=platform.platform(), vcpu=os.cpu_count())
    q = sh(['nvidia-smi', '--query-gpu=name,driver_version,memory.total,compute_cap', '--format=csv,noheader'])
    if q.returncode == 0 and q.stdout.strip():
        lines = q.stdout.strip().splitlines()
        n, d, mem, cc = [x.strip() for x in lines[0].split(',')]
        m.update(gpu=n, driver=d, gpu_mem=mem, compute_cap=cc, gpu_count=len(lines))
    nv = sh(['nvcc', '--version']) if shutil_which('nvcc') else None
    mm = re.search(r'release ([\d.]+)', nv.stdout) if nv and nv.returncode == 0 else None
    m['cuda_toolkit'] = mm.group(1) if mm else 'nvcc not found'
    smi = sh(['nvidia-smi']) if shutil_which('nvidia-smi') else None
    mm = re.search(r'CUDA Version:\s*([\d.]+)', smi.stdout) if smi else None
    m['cuda_driver_api'] = mm.group(1) if mm else 'unknown'
    try:
        m['ram_gb'] = round(int(re.search(r'MemTotal:\s+(\d+)', open('/proc/meminfo').read()).group(1)) / 1048576, 1)
    except (OSError, AttributeError):
        m['ram_gb'] = None
    cg = sh(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'])
    m['repo_commit'] = cg.stdout.strip() if cg.returncode == 0 else 'unknown'
    try:
        import numpy
        m['numpy'] = numpy.__version__
    except ImportError:
        m['numpy'] = None
    return m


def shutil_which(x):
    import shutil
    return shutil.which(x)


def build(work, mach):
    exe = os.path.join(work, 'pdhg')
    if os.path.exists(exe):
        return exe
    cc = mach.get('compute_cap', '').replace('.', '')
    arch = f'sm_{cc}' if cc.isdigit() else 'sm_75'
    cmd = ['nvcc', '-O3', f'-arch={arch}', '-std=c++17', '-Xcompiler', '-O3 -march=native -pthread', '-o', exe + '.tmp',
           os.path.join(HERE, 'pdhg.cu'), os.path.join(ROOT, 'src', 'mps.cpp')]
    print('build:', ' '.join(cmd), flush=True)
    r = sh(cmd)
    if r.returncode:
        raise SystemExit('build failed:\n' + r.stderr[-2000:])
    os.replace(exe + '.tmp', exe)
    mach['build_arch'] = arch
    return exe


def make_model(r, models):
    path = os.path.join(models, case_name(r) + '.mps')
    if os.path.exists(path):
        return path
    gen = [sys.executable, os.path.join(HERE, 'gen_lp.py')]
    gen += ['transport', *r['size'].split('x')] if r['family'] == 'transport' else ['packing', r['size']]
    t = time.time()
    p = sh(gen + [path + '.tmp'])
    if p.returncode:
        raise RuntimeError('generator failed: ' + p.stderr[-300:])
    os.replace(path + '.tmp', path)
    print(f'  generated {os.path.basename(path)} in {time.time() - t:.0f} s', flush=True)
    return path


def round2_refs():
    refs = {}
    for f in sorted(os.listdir(ROUND2)) if os.path.isdir(ROUND2) else []:
        if f.startswith('synth_') and f.endswith('.csv'):
            for row in csv.DictReader(open(os.path.join(ROUND2, f))):
                if row.get('highs_objective'):
                    refs[row['case']] = float(row['highs_objective'])
    return refs


def merge(rows_dir, cases, out):
    from bench_pdhg import FIELDS
    groups = {}
    for r in cases:
        base = os.path.join(rows_dir, r['id'])
        if os.path.exists(base + '.csv') and os.path.exists(base + '.json'):
            groups.setdefault((r['family'], r['tol']), []).append(base)
    for (fam, tol), bases in groups.items():
        stem = os.path.join(out, f'synth_{fam}_tol{tol}')
        allrows, jr, js = [], [], None
        for b in bases:
            allrows += list(csv.DictReader(open(b + '.csv')))
            js = json.load(open(b + '.json'))
            jr += js['rows']
        with open(stem + '.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
            w.writeheader()
            w.writerows(allrows)
        js['rows'] = jr
        json.dump(js, open(stem + '.json', 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', default=os.path.join(ROOT, 'out', 'round3_work'))
    ap.add_argument('--out', default=os.path.join(ROOT, 'out', 'round3'))
    ap.add_argument('--budget-s', type=float, default=6600, help='do not start a row whose estimate would pass this many seconds')
    ap.add_argument('--only', default='', help='comma list of row ids')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()
    cases = parse_cases(CASES_MD)
    if a.only:
        cases = [c for c in cases if c['id'] in a.only.split(',')]
    cases.sort(key=lambda c: c['est'])  # cheap controls first, so a budget cut-off drops the biggest rows
    if a.dry_run:
        for c in cases:
            print(c['id'], case_name(c), c['tol'], c['configs'], c['cap'], c['repeats'], c['check'], c['est'])
        print('estimated solver seconds:', sum(c['est'] for c in cases))
        return 0
    for d in ('rows', 'models'):
        os.makedirs(os.path.join(a.work, d), exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    mach = machine()
    if 'gpu' not in mach:
        raise SystemExit('no GPU visible to nvidia-smi')
    if not mach['numpy']:
        raise SystemExit('numpy is required (pip install numpy)')
    exe = build(a.work, mach)
    os.environ['TARAL_PDHG'] = exe
    refs = round2_refs()
    ref_json = os.path.join(a.work, 'refs.json')
    json.dump({case_name(c): refs.get(case_name(c)) for c in cases}, open(ref_json, 'w'))
    state_path = os.path.join(a.work, 'state.json')
    state = json.load(open(state_path)) if os.path.exists(state_path) else {}
    t_start = time.time()
    for c in cases:
        base = os.path.join(a.work, 'rows', c['id'])
        if os.path.exists(base + '.csv') and os.path.exists(base + '.json'):
            print(f"{c['id']}: finished earlier, skipped", flush=True)
            continue
        if time.time() - t_start + c['est'] > a.budget_s:
            state[c['id']] = dict(outcome='skipped: budget', est_s=c['est'])
            print(f"{c['id']}: skipped, budget ({a.budget_s:.0f} s) would be exceeded", flush=True)
            continue
        print(f"{c['id']}: {case_name(c)} tol {c['tol']} configs {c['configs']} repeats {c['repeats']} cap {c['cap']:.0f} s", flush=True)
        t0 = time.time()
        try:
            model = make_model(c, os.path.join(a.work, 'models'))
            cfg = ['gpu'] + ([f"cpu{mach['vcpu']}"] if 'cpu' in c['configs'] else [])
            check, note = c['check'], ''
            if check and os.path.getsize(model) > CHECK_MODEL_BYTES and (mach.get('ram_gb') or 0) < CHECK_RAM_GB:
                check, note = False, f'vector check skipped: model over 400 MB and RAM under {CHECK_RAM_GB} GB'
            elif not check:
                note = 'vector check not run (case list says no)'
            cmd = [sys.executable, os.path.join(HERE, 'bench_pdhg.py'), model, '--out', base + '.part', '--configs', ','.join(cfg),
                   '--tol', c['tol'], '--repeats', str(c['repeats']), '--repeat-all', '--time-limit', str(c['cap']),
                   '--ref-json', ref_json] + ([] if check else ['--no-check'])
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode or not os.path.exists(base + '.part.csv'):
                raise RuntimeError((p.stderr or p.stdout)[-400:])
            for ext in ('csv', 'json'):
                os.replace(f'{base}.part.{ext}', f'{base}.{ext}')
            rows = list(csv.DictReader(open(base + '.csv')))
            state[c['id']] = dict(outcome='finished', wall_s=round(time.time() - t0, 1), note=note,
                                  statuses={f"{r['device']}{r['threads']}": r['status'] for r in rows})
            g = [r for r in rows if r['device'] == 'gpu']
            if refs.get(case_name(c)) is not None and g and g[0]['objective_rel_error_vs_highs'] \
                    and float(g[0]['objective_rel_error_vs_highs']) > 1e-2:
                state[c['id']]['warning'] = 'objective differs from the round-2 reference by over 1e-2: model mismatch?'
            print(f"  done in {time.time() - t0:.0f} s (estimate {c['est']:.0f} s): {state[c['id']]['statuses']}", flush=True)
        except Exception as e:  # recorded, the next row still runs; a re-run retries this one
            state[c['id']] = dict(outcome='failed', error=str(e)[:400], wall_s=round(time.time() - t0, 1))
            print(f"  FAILED: {str(e)[:200]}", flush=True)
        json.dump(state, open(state_path, 'w'), indent=1)
    json.dump(state, open(state_path, 'w'), indent=1)
    merge(os.path.join(a.work, 'rows'), cases, a.out)
    write_notes(a.out, mach, cases, state, refs, time.time() - t_start, a)
    print('ledgers and NOTES.md in', a.out)
    return 0


def write_notes(out, mach, cases, state, refs, wall, a):
    L = ['# Round-3 GPU PDHG run (cloud box; measured, near-optimal method, not an exact solve)', '',
         'Case list: `gpu/ROUND3_CASES.md` (each row justified by a round-2 number). Harness: `gpu/run_round3.py`. '
         'Ledgers in the schema of `results/gpu_pdhg/synth_*.csv`: `synth_<family>_tol<tol>.csv/json`.', '', '## Environment', '']
    for k in ('gpu', 'gpu_count', 'gpu_mem', 'compute_cap', 'build_arch', 'driver', 'cuda_toolkit', 'cuda_driver_api', 'vcpu', 'ram_gb',
              'platform', 'python', 'numpy', 'repo_commit'):
        if k in mach:
            L.append(f'- {k}: {mach[k]}')
    L += [f'- harness wall time this invocation: {wall:.0f} s (budget {a.budget_s:.0f} s; resumed invocations add up)', '',
          '## Rows', '', '| id | case | tol | outcome | wall s (estimate s) | notes |', '|---|---|---|---|---|---|']
    for c in cases:
        s = state.get(c['id'], dict(outcome='not run'))
        notes = '; '.join(x for x in (s.get('note'), s.get('warning'), s.get('error')) if x)
        L.append(f"| {c['id']} | {case_name(c)} | {c['tol']} | {s['outcome']} | {s.get('wall_s', '-')} ({c['est']:.0f}) | {notes} |")
    bad = [c['id'] for c in cases if state.get(c['id'], {}).get('outcome') != 'finished']
    used = sorted(k for k, v in refs.items() if v is not None and any(case_name(c) == k for c in cases))
    L += ['', '## Failures, skips and deviations', '',
          '- Rows not finished: ' + (', '.join(bad) + ' (re-run the harness to retry; finished rows are kept)' if bad else 'none'),
          '- No CPU simplex or HiGHS was run on this box. HiGHS objectives are the round-2 values for models regenerated with the same seed ('
          + ', '.join(used) + '); other sizes have no reference and show no objective error.',
          '- Per-row detail (solver status per config, repeats, device memory, vector-check fields) is in the ledgers. '
          'A status starting with time_limit is not converged; read its error columns accordingly.',
          '- CPU PDHG rows use all vCPUs of this box (recorded as `threads`), not the 10 threads of round 2.']
    open(os.path.join(out, 'NOTES.md'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    sys.exit(main())
