"""usage: drv.py BIN LABEL tol case... ; appends rows to res_LABEL.jsonl"""
import sys, json, os, subprocess, time, numpy as np, scipy.sparse as sp, highspy
D = os.path.dirname(os.path.abspath(__file__)); BIN, LABEL, TOL = sys.argv[1:4]; cases = sys.argv[4:]
ref = json.load(open(f'{D}/ref.json'))
def load(name):
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.readModel(f'{D}/data/mps/{name}.mps')
    lp = h.getLp(); a = lp.a_matrix_
    A = sp.csc_matrix((a.value_, a.index_, a.start_), shape=(lp.num_row_, lp.num_col_))
    names = [h.getColName(j)[1] for j in range(lp.num_col_)]
    return lp, A, names
def readvec(p):
    d = {}
    for l in open(p):
        t = l.split()
        if len(t) >= 2: d[t[0]] = float(t[1])
    return d
for name in cases:
    sol = f'{D}/sol_{LABEL}_{name}'; js = sol + '.json'
    t = time.time()
    subprocess.run([BIN, f'{D}/data/mps/{name}.mps', '--device', 'gpu', '--tol', TOL, '--time-limit', '300', '--json', js, '--sol', sol], check=True, stdout=subprocess.DEVNULL)
    wall = time.time() - t; r = json.load(open(js))
    lp, A, names = load(name); v = readvec(sol)
    x = np.array([v.get(n, 0.0) for n in names]); ax = A @ x
    lo, up = np.array(lp.row_lower_), np.array(lp.row_upper_)
    viol = np.maximum(np.maximum(lo - ax, ax - up), 0)
    b = np.where(lo - ax > ax - up, lo, up); b = np.where(np.isfinite(b), b, 0)
    rel = viol / (1 + np.abs(b))
    bv = max(np.max(np.maximum(np.array(lp.col_lower_) - x, 0)), np.max(np.maximum(x - np.array(lp.col_upper_), 0)))
    obj = float(np.dot(lp.col_cost_, x)) + lp.offset_
    ho = ref[name]['objective']
    row = dict(case=name, tol=TOL, label=LABEL, status=r.get('status'), iters=r.get('iterations'), solve_s=r.get('solve_s'), setup_s=r.get('setup_s'),
               total_s=r.get('total_s'), wall_s=wall, obj=obj, ref=ho, obj_err=abs(obj - ho) / max(1, abs(ho)), worst_row=float(rel.max()), bound_viol=float(bv))
    print(json.dumps(row), flush=True); open(f'{D}/res_{LABEL}.jsonl', 'a').write(json.dumps(row) + '\n')
