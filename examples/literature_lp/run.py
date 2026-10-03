"""Measure the two published LP fixtures; retain points and independent checks.
Run: python examples/literature_lp/run.py --engine /path/to/taral --out DIR
Needs highspy and numpy for benchmark checks only, not the solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from orig_check import parse_mps, check
import highspy

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--engine', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    engine = str(Path(a.engine).resolve())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in ('steel_production', 'transportation'):
        model = Path(__file__).parent / (name + '.mps')
        sol, js = out / (name + '.sol'), out / (name + '.json')
        start = time.perf_counter()
        proc = subprocess.run([engine, str(model), '--time-limit', '30', '--sol', str(sol), '--json', str(js)], capture_output=True, text=True)
        wall = time.perf_counter() - start
        (out / (name + '.stdout.txt')).write_text(proc.stdout)
        (out / (name + '.stderr.txt')).write_text(proc.stderr)
        result = json.loads(js.read_text())
        values = dict((v[0], float(v[1])) for line in sol.read_text().splitlines() if len(v := line.split()) == 2)
        verified = check(parse_mps(model), values)
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', 1)
        assert h.readModel(str(model)) == highspy.HighsStatus.kOk
        h.run(); ref = h.getInfo().objective_function_value
        error = abs(verified['objective'] - ref) / max(1.0, abs(ref))
        passed = (proc.returncode == 0 and result['status'] == 'optimal' and h.getModelStatus() == highspy.HighsModelStatus.kOptimal
                  and verified['row_violation_rel'] <= 1e-8 and verified['bound_violation'] <= 1e-8
                  and error <= 1e-8 and abs(result['objective'] - verified['objective']) <= 1e-8 * max(1.0, abs(ref)))
        rows.append(dict(case=name, model_sha256=sha(model), status=result['status'], objective=result['objective'],
                         process_wall_seconds=wall, engine_wall_seconds=result.get('wall'), iterations=result.get('iterations'),
                         point=values, independent_original_model_check=verified, reference_status=h.modelStatusToString(h.getModelStatus()),
                         reference_objective=ref, relative_objective_error=error, passed=passed))
        print(name, result['status'], result['objective'], 'verified', passed)
        assert passed
    data = dict(source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                source_dirty=subprocess.check_output(['git', 'status', '--porcelain', 'src'], cwd=ROOT, text=True).strip(),
                src_file_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'src').glob('*')) if p.is_file()},
                engine_binary_sha256=sha(engine), highspy_version=highspy.Highs().version(), engine_cap_seconds=30, cases=rows)
    (out/'results.json').write_text(json.dumps(data, indent=2)+'\n')

if __name__ == '__main__':
    main()
