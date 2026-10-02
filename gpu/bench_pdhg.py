#!/usr/bin/env python3
"""Benchmark driver for the TARAL-LP PDHG prototype (gpu/pdhg.cu). Measured, same machine.

For each MPS case: runs the PDHG binary on every requested device config (median of --repeats runs of setup_s,
solve_s, total_s; status/objective/iterations from the first run, which is deterministic per config), solves the
ORIGINAL MPS with HiGHS (highspy, tooling only; cached in out/gpu/highs_cache.json) and writes the objective
relative error |obj - highs| / max(1, |highs|). Optionally runs the exact C++ simplex (src/) for context.
PDHG is near-optimal to a tolerance, not exact: never compare its objective with an exact-solve headline.

Usage: ~/.venvs/taral-gpu/bin/python gpu/bench_pdhg.py CASE.mps... --out out/gpu/NAME [--configs gpu,cpu1,cpu10]
       [--tol 1e-4] [--repeats 3] [--time-limit 60] [--simplex] [--fp32]
Writes NAME.csv and NAME.json (one row per case x config).
"""
import argparse, csv, json, os, statistics, subprocess, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vector_check

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDHG = os.path.join(ROOT, 'out', 'gpu', 'pdhg')
SIMPLEX = os.path.join(ROOT, 'out', 'gpu', 'taral_simplex')
CACHE = os.path.join(ROOT, 'out', 'gpu', 'highs_cache.json')
FIELDS = ['case', 'rows', 'cols', 'nnz', 'device', 'threads', 'precision', 'tol', 'status', 'objective', 'highs_objective',
          'objective_rel_error_vs_highs', 'rel_primal', 'rel_dual', 'rel_gap', 'max_row_violation', 'iterations',
          'restarts', 'setup_s', 'solve_s', 'total_s', 'total_s_min', 'total_s_max', 'repeats', 'device_mem_bytes',
          'highs_s', 'load_avg_1m', 'chk_objective', 'chk_row_viol_rel', 'chk_bound_viol', 'chk_dual_bound', 'chk_dual_obj_loose',
          'chk_dual_infeas_rel', 'chk_dual_gap_rel', 'chk_dual_sign', 'chk_row_viol_rel_max_repeats', 'objective_identical_repeats', 'chk_error']


def highs_ref(path, cache, solver=None):
    key = os.path.abspath(path)
    st = os.stat(path)
    hit = cache.get(key)
    if hit and hit.get('mtime') == st.st_mtime and hit.get('size') == st.st_size:
        return hit
    import highspy
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    if solver:
        h.setOptionValue('solver', solver)
    h.readModel(path)
    t = time.perf_counter()
    h.run()
    wall = time.perf_counter() - t
    status = h.modelStatusToString(h.getModelStatus())
    obj = h.getInfo().objective_function_value if status == 'Optimal' else None
    cache[key] = hit = dict(status=status, objective=obj, wall_s=wall, solver=solver or 'choose', mtime=st.st_mtime, size=st.st_size)
    if os.path.exists(CACHE):  # merge: another bench process may have added entries meanwhile
        cache.update({k: v for k, v in json.load(open(CACHE)).items() if k not in cache})
    with open(CACHE + '.tmp', 'w') as f:
        json.dump(cache, f, indent=1)
    os.replace(CACHE + '.tmp', CACHE)
    return hit


def wait_load(a):
    t = time.time()
    while a.max_load > 0 and os.getloadavg()[0] > a.max_load and time.time() - t < 1200:
        time.sleep(10)


def run_pdhg(path, device, threads, a):
    wait_load(a)
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as tf:
        out = tf.name
    sol, dual = out + '.sol', out + '.dual'
    cmd = [PDHG, path, '--device', device, '--tol', str(a.tol), '--time-limit', str(a.time_limit), '--json', out,
           '--sol', sol, '--dual', dual]
    if threads:
        cmd += ['--threads', str(threads)]
    if a.fp32:
        cmd += ['--fp32']
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=a.time_limit + 120, check=True)
        with open(out) as f:
            r = json.load(f)
        # Independent recheck of the written vectors on the ORIGINAL model (gpu/vector_check.py).
        try:
            r.update(vector_check.check_files(path, sol, dual if os.path.exists(dual) else None))
        except Exception as e:  # a failed check is recorded, never silently skipped
            r['chk_error'] = str(e)[:200]
        return r
    except (subprocess.SubprocessError, OSError, json.JSONDecodeError) as e:
        return dict(status='crash', error=str(e)[:200])
    finally:
        for f in (out, sol, dual):
            if os.path.exists(f):
                os.unlink(f)


def run_simplex(path, a):
    wait_load(a)
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as tf:
        out = tf.name
    try:
        t = time.perf_counter()
        subprocess.run([SIMPLEX, path, '--time-limit', str(a.time_limit), '--json', out], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=a.time_limit + 120)
        wall = time.perf_counter() - t
        with open(out) as f:
            d = json.load(f)
        d['process_wall_s'] = wall
        return d
    except (subprocess.SubprocessError, OSError, json.JSONDecodeError) as e:
        return dict(status='crash', error=str(e)[:200])
    finally:
        if os.path.exists(out):
            os.unlink(out)


def rel_err(obj, ref):
    if obj is None or ref is None:
        return None
    return abs(obj - ref) / max(1.0, abs(ref))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cases', nargs='+')
    ap.add_argument('--out', required=True)
    ap.add_argument('--configs', default='gpu,cpu1,cpu10', help='gpu and/or cpuN (N threads)')
    ap.add_argument('--tol', type=float, default=1e-4)
    ap.add_argument('--repeats', type=int, default=3)
    ap.add_argument('--time-limit', type=float, default=60)
    ap.add_argument('--simplex', action='store_true', help='also run the exact C++ simplex for context')
    ap.add_argument('--fp32', action='store_true')
    ap.add_argument('--repeat-all', action='store_true', help='repeat every non-crashed config, also time-limited ones')
    ap.add_argument('--max-load', type=float, default=0, help='wait (up to 20 min) for 1-min load below this')
    a = ap.parse_args()
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    rows = []
    for path in a.cases:
        case = os.path.splitext(os.path.basename(path))[0]
        ref = highs_ref(path, cache)
        for cfg in a.configs.split(','):
            device, threads = ('gpu', 0) if cfg == 'gpu' else ('cpu', int(cfg[3:]))
            runs = [run_pdhg(path, device, threads, a)]
            # Repeat only finished runs: a time-limited or crashed config is not re-timed.
            if runs[0].get('status') == 'near_optimal' or (a.repeat_all and runs[0].get('status') != 'crash'):
                runs += [run_pdhg(path, device, threads, a) for _ in range(a.repeats - 1)]
            r0 = runs[0]
            ok = [r for r in runs if 'total_s' in r]
            med = lambda k: statistics.median(r[k] for r in ok) if ok else None
            row = dict(case=case, rows=r0.get('rows'), cols=r0.get('cols'), nnz=r0.get('nnz'), device=device,
                       threads=r0.get('threads'), precision=r0.get('precision'), tol=a.tol, status=r0.get('status'),
                       objective=r0.get('objective'), highs_objective=ref['objective'],
                       objective_rel_error_vs_highs=rel_err(r0.get('objective'), ref['objective']),
                       rel_primal=r0.get('rel_primal'), rel_dual=r0.get('rel_dual'), rel_gap=r0.get('rel_gap'),
                       max_row_violation=r0.get('max_row_violation'), iterations=r0.get('iterations'),
                       restarts=r0.get('restarts'), setup_s=med('setup_s'), solve_s=med('solve_s'),
                       total_s=med('total_s'), total_s_min=min((r['total_s'] for r in ok), default=None),
                       total_s_max=max((r['total_s'] for r in ok), default=None), repeats=len(ok),
                       device_mem_bytes=r0.get('device_mem_bytes'), highs_s=ref['wall_s'],
                       load_avg_1m=os.getloadavg()[0],
                       **{k: r0.get(k) for k in ('chk_objective', 'chk_row_viol_rel', 'chk_bound_viol', 'chk_dual_bound',
                                                'chk_dual_obj_loose', 'chk_dual_infeas_rel', 'chk_dual_gap_rel', 'chk_dual_sign')},
                       chk_row_viol_rel_max_repeats=max((r.get('chk_row_viol_rel') or 0 for r in ok), default=None),
                       objective_identical_repeats=len({r.get('objective') for r in ok}) <= 1 if ok else None)
            errs = sorted({r['chk_error'] for r in runs if r.get('chk_error')})
            if errs:  # a run whose vectors could not be rechecked is vetoed, never reported as its solver status
                row['chk_error'] = '; '.join(errs)
                row['status'] = 'chk_failed (solver said %s)' % row['status']
            if len({r.get('iterations') for r in ok}) > 1:
                row['status'] += ' (iterations differ across repeats)'
            rows.append(row)
            e = row['objective_rel_error_vs_highs']
            print(f"{case:12s} {cfg:6s} {row['status']:14s} it {row['iterations']} relerr "
                  f"{'n/a' if e is None else f'{e:.2e}'} total {row['total_s']} solve {row['solve_s']}", flush=True)
        if a.simplex:
            s = run_simplex(path, a)
            obj = s.get('objective')
            rows.append(dict(case=case, rows=rows[-1]['rows'], cols=rows[-1]['cols'], nnz=rows[-1]['nnz'],
                             device='cpu-exact-simplex', threads=1, precision='fp64', tol='exact',
                             status=s.get('status'), objective=obj, highs_objective=ref['objective'],
                             objective_rel_error_vs_highs=rel_err(obj, ref['objective']),
                             iterations=s.get('iterations'), total_s=s.get('wall_s'), repeats=1,
                             highs_s=ref['wall_s'], load_avg_1m=os.getloadavg()[0]))
            print(f"{case:12s} simplex {s.get('status')} total {s.get('wall_s')}", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out + '.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    with open(a.out + '.json', 'w') as f:
        json.dump(dict(tool='taral-pdhg bench (prototype, near-optimal; measured)', tol=a.tol, repeats=a.repeats,
                       time_limit=a.time_limit, fp32=a.fp32, rows=rows), f, indent=1)


if __name__ == '__main__':
    sys.exit(main())
