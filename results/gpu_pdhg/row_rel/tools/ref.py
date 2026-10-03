import sys, json, time, os, highspy
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ref.json')
ref = json.load(open(out)) if os.path.exists(out) else {}
for name in sys.argv[1:]:
    if name in ref: continue
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'ipm'); h.setOptionValue('threads', 4)
    h.setOptionValue('time_limit', 3000.0)
    h.readModel(f'data/mps/{name}.mps'); t = time.time(); h.run()
    st = h.modelStatusToString(h.getModelStatus()); ob = h.getInfo().objective_function_value
    ref[name] = dict(status=st, objective=ob, s=time.time()-t); print(name, ref[name], flush=True)
    json.dump(ref, open(out, 'w'), indent=1)
