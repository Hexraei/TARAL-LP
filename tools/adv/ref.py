"""HiGHS reference (highspy) for the adversarial suite. HiGHS is used here only; the engine never links it."""
import numpy as np
import highspy

from mpsio import INF

STATUS = {
    highspy.HighsModelStatus.kOptimal: "optimal",
    highspy.HighsModelStatus.kInfeasible: "infeasible",
    highspy.HighsModelStatus.kUnbounded: "unbounded",
    highspy.HighsModelStatus.kUnboundedOrInfeasible: "unbounded_or_infeasible",
    highspy.HighsModelStatus.kTimeLimit: "time_limit",
    highspy.HighsModelStatus.kModelEmpty: "optimal",  # empty model: trivially optimal at the origin
}


def _lp(m):
    n, r = len(m.col_names), len(m.row_names)
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = n, r
    lp.col_cost_ = np.array(m.cost, dtype=float)
    lp.col_lower_ = np.array(m.col_lo, dtype=float)
    lp.col_upper_ = np.array(m.col_up, dtype=float)
    lp.row_lower_ = np.array(m.row_lo, dtype=float)
    lp.row_upper_ = np.array(m.row_up, dtype=float)
    start, index, value = [0], [], []
    for j in range(n):
        for i in sorted(m.entries[j]):
            v = m.entries[j][i]
            if v != 0.0:
                index.append(i)
                value.append(v)
        start.append(len(index))
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = np.array(start, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(index, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(value, dtype=float)
    lp.offset_ = m.obj_const
    lp.sense_ = highspy.ObjSense.kMaximize if m.maximize else highspy.ObjSense.kMinimize
    if any(m.is_int):
        lp.integrality_ = [highspy.HighsVarType.kInteger if v else highspy.HighsVarType.kContinuous for v in m.is_int]
    return lp


def _solve_ref(m, presolve=True, time_limit=60.0, solver="choose"):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("presolve", "on" if presolve else "off")
    h.setOptionValue("solver", solver)
    h.setOptionValue("threads", 1)
    h.setOptionValue("random_seed", 0)
    # feasibility tolerances stay at the HiGHS defaults: tightening mip_feasibility_tolerance to 1e-9 made HiGHS 1.15.1
    # cut off the true optimum of a 1e7-scaled MIP (see tests/repro/big_values_numfail.mps)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 0.0)
    h.setOptionValue("small_matrix_value", 1e-12)  # default 1e-9 silently drops 1e-9-sized coefficients
    h.passModel(_lp(m))
    h.run()
    ms = h.getModelStatus()
    st = STATUS.get(ms, "other:" + str(ms).split(".")[-1])
    out = dict(status=st, objective=None, x=None)
    if st == "optimal":
        out["objective"] = h.getInfo().objective_function_value
        out["x"] = list(h.getSolution().col_value)
        if not out["x"]:  # model with zero columns
            out["x"] = []
    return out


def read_ref(path):
    """What HiGHS makes of the file itself (used by the parser table and the load-agreement check)."""
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    st = h.readModel(path)
    ok = st != highspy.HighsStatus.kError
    return h, ok, str(st).split(".")[-1]


def solve_ref(m, presolve=True, time_limit=60.0, solver="choose", hard=None):
    """_solve_ref in a forked child with a hard wall-clock cap: HiGHS 1.15.1 can ignore time_limit (hang) or crash on
    some badly scaled models, and neither may take the harness down. Status 'ref_hang' / 'ref_crash' is reported."""
    import multiprocessing as mp
    hard = hard if hard is not None else time_limit + 15
    ctx = mp.get_context("fork")
    parent, child = ctx.Pipe(duplex=False)

    def run():
        try:
            child.send(_solve_ref(m, presolve, time_limit, solver))
        except Exception as e:  # pragma: no cover
            child.send(dict(status="ref_exception:" + repr(e)[:60], objective=None, x=None))

    p = ctx.Process(target=run)
    p.start()
    child.close()
    if parent.poll(hard):
        try:
            out = parent.recv()
        except EOFError:
            out = dict(status="ref_crash", objective=None, x=None)
    else:
        out = dict(status="ref_hang", objective=None, x=None)
    if p.is_alive():
        p.kill()
    p.join()
    if out is None:
        out = dict(status="ref_crash", objective=None, x=None)
    return out
