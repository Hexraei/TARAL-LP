#!/usr/bin/env python3
"""Delta-debugging shrinker: reduce a failing generated case to a small reproducer.

The failure signature is (engine status, reference status) with the reference solved by HiGHS with presolve off,
or 'objective' when both are optimal but the objectives differ by more than 1e-6 relative. The reduced model is
re-written as clean free-format MPS (R<i>/C<i> names) from the oracle's parsed model, so the reproducer does not
depend on any parser quirk of the original file.

Usage: shrink.py --engine ./taral --case cat-seed --out tests/repro/NAME.mps [--time-limit 3] [--method simplex|dual|ipm]
"""
import argparse
import copy
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from mpsio import parse_mps, num, INF, Model
from ref import solve_ref


def write_model(m, path, sentinel=0.0):
    """Clean free-format MPS from an oracle Model. Rows with both sides infinite are dropped.
    sentinel=S (e.g. 1e20 or 1e30) writes every infinite column bound as +-S (so sentinel-handling bugs survive shrinking)."""
    rows = [i for i in range(len(m.row_names)) if not (m.row_lo[i] == -INF and m.row_up[i] == INF)]
    out = ["NAME          REPRO"]
    if m.maximize:
        out += ["OBJSENSE", "    MAX"]
    out.append("ROWS")
    out.append(" N  OBJ")
    ranges, rhs = [], []
    for i in rows:
        lo, up = m.row_lo[i], m.row_up[i]
        nm = "R%d" % i
        if lo == up:
            out.append(" E  " + nm), rhs.append((nm, lo))
        elif lo == -INF:
            out.append(" L  " + nm), rhs.append((nm, up))
            if sentinel:
                ranges.append((nm, sentinel))  # one-sided row spelled with an infinite-sentinel range
        elif up == INF:
            out.append(" G  " + nm), rhs.append((nm, lo))
            if sentinel:
                ranges.append((nm, sentinel))
        else:
            out.append(" G  " + nm), rhs.append((nm, lo)), ranges.append((nm, up - lo))
    out.append("COLUMNS")
    in_int = False
    mk = 0
    for j in range(len(m.col_names)):
        if bool(m.is_int[j]) != in_int:
            out.append("    MARKER%d  'MARKER'  '%s'" % (mk, "INTORG" if m.is_int[j] else "INTEND"))
            mk += 1
            in_int = bool(m.is_int[j])
        ent = []
        if m.cost[j] != 0:
            ent.append(("OBJ", m.cost[j]))
        ent += [("R%d" % i, v) for i, v in sorted(m.entries[j].items()) if i in set(rows) and v != 0]
        if not ent:
            ent = [("OBJ", 0)]
        for k in range(0, len(ent), 2):
            out.append("    C%d  " % j + "  ".join("%s  %s" % (r, num(v)) for r, v in ent[k:k + 2]))
    if in_int:
        out.append("    MARKER%d  'MARKER'  'INTEND'" % mk)
    if m.obj_const != 0:
        rhs.insert(0, ("OBJ", -m.obj_const))
    if rhs:
        out.append("RHS")
        out += ["    RHS  %s  %s" % (r, num(v)) for r, v in rhs]
    if ranges:
        out.append("RANGES")
        out += ["    RNG  %s  %s" % (r, num(v)) for r, v in ranges]
    bd = []
    for j in range(len(m.col_names)):
        lo, up, c = m.col_lo[j], m.col_up[j], "C%d" % j
        if sentinel and (lo == -INF or up == INF):
            lo = -sentinel if lo == -INF else lo
            up = sentinel if up == INF else up
            if lo != 0:
                bd.append(" LO BND  %s  %s" % (c, num(lo)))
            bd.append(" UP BND  %s  %s" % (c, num(up)))
            continue
        if lo == 0 and up == INF:
            if m.is_int[j]:
                bd.append(" PL BND  " + c)
            continue
        if lo == -INF and up == INF:
            bd.append(" FR BND  " + c)
        elif lo == up:
            bd.append(" FX BND  %s  %s" % (c, num(lo)))
        else:
            if lo == -INF:
                bd.append(" MI BND  " + c)
            elif lo != 0:
                bd.append(" LO BND  %s  %s" % (c, num(lo)))
            if up != INF:
                bd.append(" UP BND  %s  %s" % (c, num(up)))
    if bd:
        out.append("BOUNDS")
        out += bd
    out.append("ENDATA")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def drop_row(m, i):
    n = copy.deepcopy(m)
    del n.row_names[i], n.row_lo[i], n.row_up[i]
    n.entries = [{(k if k < i else k - 1): v for k, v in e.items() if k != i} for e in n.entries]
    return n


def drop_col(m, j):
    n = copy.deepcopy(m)
    for lst in (n.col_names, n.cost, n.col_lo, n.col_up, n.is_int, n.entries):
        del lst[j]
    return n


class Shrinker:
    def __init__(self, engine, tl, mip, sentinel=0.0, method="simplex"):
        self.engine, self.tl, self.mip, self.sentinel, self.method = engine, tl, mip, sentinel, method
        self.tmp = tempfile.mkdtemp(prefix="shrink")
        self.calls = 0

    def signature(self, m):
        self.calls += 1
        p = os.path.join(self.tmp, "c.mps")
        js = p + ".json"
        write_model(m, p, self.sentinel)
        if os.path.exists(js):
            os.remove(js)
        try:
            cmd = [self.engine, p, "--time-limit", str(self.tl), "--json", js, "--sol", p + ".sol"]
            if self.method != "simplex":
                cmd += ["--method", self.method]
            subprocess.run(cmd, capture_output=True, timeout=self.tl + 20)
            j = json.load(open(js))
        except Exception:
            return None
        ours = j.get("status")
        ref = solve_ref(parse_mps(p), presolve=False, time_limit=30)
        if ref["status"] not in ("optimal", "infeasible", "unbounded"):
            return None
        if ours != ref["status"]:
            return ("status", ours, ref["status"])
        if ours == "optimal" and j.get("objective") is not None and ref["objective"] is not None:
            if abs(j["objective"] - ref["objective"]) > 1e-6 * max(1.0, abs(ref["objective"])):
                return ("objective",)
        return None

    def run(self, m, sig):
        ok = lambda x: self.signature(x) == sig
        changed = True
        while changed:
            changed = False
            for i in reversed(range(len(m.row_names))):
                c = drop_row(m, i)
                if ok(c):
                    m, changed = c, True
            for j in reversed(range(len(m.col_names))):
                if len(m.col_names) > 1:
                    c = drop_col(m, j)
                    if ok(c):
                        m, changed = c, True
            for j in range(len(m.col_names)):
                for i in sorted(m.entries[j]):
                    c = copy.deepcopy(m)
                    del c.entries[j][i]
                    if ok(c):
                        m, changed = c, True
                if m.cost[j] != 0:
                    c = copy.deepcopy(m)
                    c.cost[j] = 0.0
                    if ok(c):
                        m, changed = c, True
            # simplify bounds and values
            for j in range(len(m.col_names)):
                for field, dflt in (("col_lo", 0.0), ("col_up", INF)):
                    if getattr(m, field)[j] != dflt:
                        c = copy.deepcopy(m)
                        getattr(c, field)[j] = dflt
                        if ok(c):
                            m, changed = c, True
                if m.is_int[j] and self.mip:
                    pass
            for i in range(len(m.row_names)):
                for field in ("row_lo", "row_up"):
                    v = getattr(m, field)[i]
                    if v not in (INF, -INF) and v != 0:
                        c = copy.deepcopy(m)
                        getattr(c, field)[i] = float(round(v)) if abs(v) < 1e6 else v
                        if getattr(c, field)[i] != v and c.row_lo[i] <= c.row_up[i] and ok(c):
                            m, changed = c, True
            if m.obj_const != 0:
                c = copy.deepcopy(m)
                c.obj_const = 0.0
                if ok(c):
                    m, changed = c, True
        return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--time-limit", type=float, default=3.0)
    ap.add_argument("--sentinel", type=float, default=0.0, help="write infinite column bounds as +-S (1e20 or 1e30)")
    ap.add_argument("--method", default="simplex", choices=["simplex", "dual", "ipm"], help="engine path to reproduce on")
    a = ap.parse_args()
    cat, seed = a.case.rsplit("-", 1)
    case = gen.gen(cat, int(seed) - gen.BASE[cat])
    tmp = tempfile.mkdtemp(prefix="orig")
    p = os.path.join(tmp, "orig.mps")
    open(p, "w").write(case.text)
    m = parse_mps(p)
    sh = Shrinker(os.path.abspath(a.engine), a.time_limit, case.mip, a.sentinel, a.method)
    sig = sh.signature(m)
    if sig is None:
        print("case does not reproduce with --time-limit %s" % a.time_limit)
        sys.exit(1)
    print("signature", sig, "start size %dx%d" % (len(m.row_names), len(m.col_names)))
    m = sh.run(m, sig)
    write_model(m, a.out, a.sentinel)
    print("reduced to %dx%d in %d engine calls -> %s" % (len(m.row_names), len(m.col_names), sh.calls, a.out))
    ref = solve_ref(m, presolve=False)
    print(json.dumps(dict(signature=sig, ref=ref["status"], ref_obj=ref["objective"])))


if __name__ == "__main__":
    main()
