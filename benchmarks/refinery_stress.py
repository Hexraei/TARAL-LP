#!/usr/bin/env python3
"""Synthetic refinery stress suite: multi-period / supply / quality LP and MILP cases.

ALL DATA HERE IS SYNTHETIC. The cases are generated from a seeded recipe and are not operational data
from any refinery; the cost and yield numbers are invented to give a plausible-looking shape. The
structure is loosely modelled on the planning-model families in examples/refinery_williams.py.

Everything is LINEAR. Blend quality rows (octane, sulfur) use FIXED component qualities, so they are
linear in volumes. Nonlinear pooling (intermediate tanks whose quality is itself a variable, giving
bilinear quality x volume terms) is NOT modelled and nothing here should be read as solving it.

Modes
  --write DIR    regenerate the MPS files (deterministic; committed copies must match)
  (default)      solve every committed case with the engine, then independently check it:
                   1. own free-format MPS parser (does not use the engine's reader)
                   2. rows, bounds, integrality and objective re-evaluated from the engine's .sol
                   3. HiGHS (scipy.optimize.linprog / milp) used ONLY as reference oracle for status/objective
                   4. a scenario report (profit, crude bought, throughput, product sales) is computed from the
                      verified point and written to --report (JSON) and printed.
Usage: python3 benchmarks/refinery_stress.py [--engine ./taral] [--dir benchmarks/refinery_stress] [--report out.json]
Exit code 1 on any WRONG case. Requires numpy and scipy (>=1.9).
"""
import argparse, json, math, os, subprocess, sys, tempfile
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix

CRUDES = {  # name: (cost, [naphtha, distillate, residue yield])  (synthetic)
    "SWEET": (78.0, (0.28, 0.42, 0.28)),
    "SOUR": (68.0, (0.20, 0.38, 0.40)),
}
OCT = {"N": 72.0, "RF": 100.0, "BU": 95.0}
SULF = {"D": 0.50, "HD": 0.05}
PRICE = {"GAS": 105.0, "DSL": 98.0, "FO": 62.0}
BU_COST, HD_COST, REFORM_YIELD, REFORM_COST = 85.0, 92.0, 0.85, 4.0
OCT_SPEC, SULF_SPEC = 87.0, 0.30
INV_COST = 0.15  # per barrel held per period, all tanks


def scenarios():
    """name -> dict. T periods; avail scale per crude; demand scale; mip flag; expected truth."""
    return {
        "base_lp_t6": dict(T=6, seed=1),
        "base_lp_t12": dict(T=12, seed=2),
        "base_lp_t24": dict(T=24, seed=3),
        "supply_cut_lp": dict(T=12, seed=2, avail={"SOUR": 0.3}),
        "sour_outage_lp": dict(T=12, seed=2, avail={"SOUR": 0.0}),
        "demand_surge_lp": dict(T=12, seed=2, demand=1.35),
        "quality_tight_lp": dict(T=12, seed=2, oct_spec=93.0, sulf_spec=0.12),
        "infeasible_supply_lp": dict(T=6, seed=1, avail={"SWEET": 0.05, "SOUR": 0.05}, demand=1.0, mindemand=True,
                                     truth="infeasible"),
        "campaign_milp_t6": dict(T=6, seed=1, mip=True),
        "campaign_milp_t8": dict(T=8, seed=4, mip=True),
        "campaign_supply_cut_milp": dict(T=8, seed=4, mip=True, avail={"SOUR": 0.3}),
        "cargo_milp_t6": dict(T=6, seed=5, mip=True, cargo=True),
    }


def num(v):
    return repr(float(v))


def build(sc):
    """Return a column-oriented model dict: cols (name -> (cost,lo,up,int)), rows (name -> (sense,rhs)), coef."""
    T, rng = sc["T"], np.random.RandomState(sc["seed"])
    mip, cargo = sc.get("mip", False), sc.get("cargo", False)
    cols, rows, coef = {}, {}, {}  # coef[(row, col)] = v

    def col(n, cost=0.0, lo=0.0, up=math.inf, integer=False):
        cols[n] = (cost, lo, up, integer)

    def row(n, sense, rhs=0.0):
        rows[n] = (sense, rhs)

    def a(r, c, v):
        coef[(r, c)] = coef.get((r, c), 0.0) + v

    base_av = {"SWEET": 14000.0, "SOUR": 20000.0}
    base_dem = {"GAS": 9000.0, "DSL": 11000.0, "FO": 5500.0}
    cap_run, tank = 30000.0, 40000.0
    cap_ref, cap_hd, cap_bu = 5000.0, 4000.0, 3000.0
    scale = sc.get("demand", 1.0)
    oct_spec, sulf_spec = sc.get("oct_spec", OCT_SPEC), sc.get("sulf_spec", SULF_SPEC)
    for t in range(T):
        wave = 1.0 + 0.15 * math.sin(2 * math.pi * (t + rng.rand()) / 6)
        for cn, (cost, y) in CRUDES.items():
            av = base_av[cn] * sc.get("avail", {}).get(cn, 1.0) * (0.85 + 0.3 * rng.rand())
            if cargo:  # integer 2,000-bbl parcels; purchase = 2000 * k
                col(f"K_{cn}_{t}", cost * 2000.0, 0, math.floor(av / 2000.0), True)
                col(f"BUY_{cn}_{t}", 0.0)
                row(f"CARGO_{cn}_{t}", "E")
                a(f"CARGO_{cn}_{t}", f"BUY_{cn}_{t}", 1.0)
                a(f"CARGO_{cn}_{t}", f"K_{cn}_{t}", -2000.0)
            else:
                col(f"BUY_{cn}_{t}", cost * (0.95 + 0.1 * rng.rand()), 0, av)
            col(f"CINV_{cn}_{t}", INV_COST, 0, tank)
            col(f"RUN_{cn}_{t}")
            r = f"CBAL_{cn}_{t}"  # opening + buy - run - closing = 0
            row(r, "E", -0.0)
            a(r, f"BUY_{cn}_{t}", 1.0); a(r, f"RUN_{cn}_{t}", -1.0); a(r, f"CINV_{cn}_{t}", -1.0)
            if t > 0:
                a(r, f"CINV_{cn}_{t-1}", 1.0)
            else:
                rows[r] = ("E", -2000.0)  # opening stock 2000 bbl
        row(f"CAP_{t}", "L", cap_run)
        for cn in CRUDES:
            a(f"CAP_{t}", f"RUN_{cn}_{t}", 1.0)
        # intermediates: N, D, R
        for k, s in enumerate(("N", "D", "R")):
            col(f"IINV_{s}_{t}", INV_COST, 0, 15000.0)
            r = f"IBAL_{s}_{t}"
            row(r, "E", 0.0)
            for cn, (_, y) in CRUDES.items():
                a(r, f"RUN_{cn}_{t}", y[k])
            a(r, f"IINV_{s}_{t}", -1.0)
            if t > 0:
                a(r, f"IINV_{s}_{t-1}", 1.0)
        # N: to reformer or to gasoline blend; D: to diesel blend; R: to fuel oil
        col(f"N_RF_{t}", REFORM_COST, 0, cap_ref)
        col(f"N_G_{t}"); col(f"RF_G_{t}"); col(f"BU_G_{t}", BU_COST, 0, cap_bu)
        col(f"D_D_{t}"); col(f"HD_D_{t}", HD_COST, 0, cap_hd)
        col(f"R_FO_{t}")
        # IBAL: production - inventory change - uses = 0 -> add uses
        a(f"IBAL_N_{t}", f"N_RF_{t}", -1.0); a(f"IBAL_N_{t}", f"N_G_{t}", -1.0)
        a(f"IBAL_D_{t}", f"D_D_{t}", -1.0)
        a(f"IBAL_R_{t}", f"R_FO_{t}", -1.0)
        row(f"RFBAL_{t}", "E"); a(f"RFBAL_{t}", f"N_RF_{t}", REFORM_YIELD); a(f"RFBAL_{t}", f"RF_G_{t}", -1.0)
        # gasoline blend, quality (fixed component qualities -> linear)
        for p, comps in (("GAS", ("N_G", "RF_G", "BU_G")), ("DSL", ("D_D", "HD_D")), ("FO", ("R_FO",))):
            col(f"SELL_{p}_{t}", -PRICE[p] * (0.97 + 0.06 * rng.rand()), 0, math.inf)
            col(f"PINV_{p}_{t}", INV_COST, 0, 12000.0)
            r = f"PBAL_{p}_{t}"
            row(r, "E", 0.0)
            for c in comps:
                a(r, f"{c}_{t}", 1.0)
            a(r, f"SELL_{p}_{t}", -1.0); a(r, f"PINV_{p}_{t}", -1.0)
            if t > 0:
                a(r, f"PINV_{p}_{t-1}", 1.0)
            dem = base_dem[p] * scale * wave
            cols[f"SELL_{p}_{t}"] = (cols[f"SELL_{p}_{t}"][0], 0.0, dem, False)
            if sc.get("mindemand"):  # contractual minimum equal to the demand: makes thin supply infeasible
                cols[f"SELL_{p}_{t}"] = (cols[f"SELL_{p}_{t}"][0], dem, dem, False)
        row(f"OCT_{t}", "G"); row(f"SULF_{t}", "L")
        for c, q in (("N_G", OCT["N"]), ("RF_G", OCT["RF"]), ("BU_G", OCT["BU"])):
            a(f"OCT_{t}", f"{c}_{t}", q - oct_spec)
        for c, q in (("D_D", SULF["D"]), ("HD_D", SULF["HD"])):
            a(f"SULF_{t}", f"{c}_{t}", q - sulf_spec)
        if mip:  # crude unit on/off with minimum stable throughput and a startup charge; reformer campaign
            col(f"ON_{t}", 12000.0, 0, 1, True)
            col(f"UP_{t}", 40000.0, 0, 1, True)
            row(f"MAXON_{t}", "L", 0.0); row(f"MINON_{t}", "G", 0.0)
            for cn in CRUDES:
                a(f"MAXON_{t}", f"RUN_{cn}_{t}", 1.0); a(f"MINON_{t}", f"RUN_{cn}_{t}", 1.0)
            a(f"MAXON_{t}", f"ON_{t}", -cap_run); a(f"MINON_{t}", f"ON_{t}", -8000.0)
            row(f"STARTUP_{t}", "L", 0.0 if t > 0 else 0.0)  # ON_t - ON_{t-1} - UP_t <= 0
            a(f"STARTUP_{t}", f"ON_{t}", 1.0); a(f"STARTUP_{t}", f"UP_{t}", -1.0)
            if t > 0:
                a(f"STARTUP_{t}", f"ON_{t-1}", -1.0)
            col(f"RFON_{t}", 5000.0, 0, 1, True)
            row(f"RFMAX_{t}", "L"); row(f"RFMIN_{t}", "G")
            a(f"RFMAX_{t}", f"N_RF_{t}", 1.0); a(f"RFMAX_{t}", f"RFON_{t}", -cap_ref)
            a(f"RFMIN_{t}", f"N_RF_{t}", 1.0); a(f"RFMIN_{t}", f"RFON_{t}", -1500.0)
    # closing stocks must at least match opening stock
    for cn in CRUDES:
        lo, up = cols[f"CINV_{cn}_{T-1}"][1:3]
        cols[f"CINV_{cn}_{T-1}"] = (cols[f"CINV_{cn}_{T-1}"][0], 2000.0, up, False)
    return cols, rows, coef


def write_mps(name, model, path):
    cols, rows, coef = model
    out = ["NAME          " + name, "* SYNTHETIC refinery stress case; invented data, linear model only.", "ROWS", " N  COST"]
    for r, (s, _) in rows.items():
        out.append(f" {s}  {r}")
    byc = {}
    for (r, c), v in coef.items():
        byc.setdefault(c, []).append((r, v))
    out.append("COLUMNS")
    inint = False
    for c, (cost, lo, up, isint) in cols.items():
        if isint and not inint:
            out.append("    MARKER                 'MARKER'                 'INTORG'"); inint = True
        if not isint and inint:
            out.append("    MARKER                 'MARKER'                 'INTEND'"); inint = False
        ent = ([("COST", cost)] if cost != 0 else []) + sorted(byc.get(c, []))
        for r, v in ent:
            out.append(f"    {c}  {r}  {num(v)}")
        if not ent:
            out.append(f"    {c}  COST  0")
    if inint:
        out.append("    MARKER                 'MARKER'                 'INTEND'")
    out.append("RHS")
    for r, (s, rhs) in rows.items():
        if rhs != 0:
            out.append(f"    RHS  {r}  {num(rhs)}")
    out.append("BOUNDS")
    for c, (cost, lo, up, isint) in cols.items():
        if lo == up:
            out.append(f" FX BND  {c}  {num(lo)}")
            continue
        if lo != 0:
            out.append(f" LO BND  {c}  {num(lo)}")
        if up != math.inf:
            out.append(f" UP BND  {c}  {num(up)}")
    out.append("ENDATA")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


# ---------------------------------------------------------------- independent checker
def parse_mps(path):
    """Minimal free-format reader for the subset this generator writes (N/L/G/E rows, MARKER, RHS, LO/UP/FX)."""
    sec, rows, rsense, cols, ints, cost = None, [], {}, [], set(), {}
    entries, rhs, lo, up, inint = [], {}, {}, {}, False
    for line in open(path):
        if line.startswith("*") or not line.strip():
            continue
        if not line.startswith(" "):
            sec = line.split()[0]
            continue
        t = line.split()
        if sec == "ROWS":
            if t[0] == "N":
                obj = t[1]
            else:
                rows.append(t[1]); rsense[t[1]] = t[0]
        elif sec == "COLUMNS":
            if len(t) >= 3 and t[1] == "'MARKER'":
                inint = t[2] == "'INTORG'"
                continue
            c = t[0]
            if c not in lo:
                cols.append(c); lo[c] = 0.0; up[c] = math.inf
                if inint:
                    ints.add(c)
            for i in range(1, len(t), 2):
                if t[i] == obj:
                    cost[c] = float(t[i + 1])
                else:
                    entries.append((t[i], c, float(t[i + 1])))
        elif sec == "RHS":
            rhs[t[1]] = float(t[2])
        elif sec == "BOUNDS":
            ty, c, v = t[0], t[2], float(t[3])
            if ty == "LO": lo[c] = v
            elif ty == "UP": up[c] = v
            elif ty == "FX": lo[c] = up[c] = v
            else: raise ValueError("unsupported bound " + ty)
    ri = {r: i for i, r in enumerate(rows)}; ci = {c: j for j, c in enumerate(cols)}
    A = np.zeros((len(rows), len(cols)))
    for r, c, v in entries:
        A[ri[r], ci[c]] += v
    b = np.array([rhs.get(r, 0.0) for r in rows])
    rlo = np.array([b[i] if rsense[r] in "GE" else -math.inf for i, r in enumerate(rows)])
    rup = np.array([b[i] if rsense[r] in "LE" else math.inf for i, r in enumerate(rows)])
    return dict(rows=rows, cols=cols, A=A, rlo=rlo, rup=rup, c=np.array([cost.get(c, 0.0) for c in cols]),
                lo=np.array([lo[c] for c in cols]), up=np.array([up[c] for c in cols]),
                isint=np.array([c in ints for c in cols]))


def highs(m):
    A = csr_matrix(m["A"])
    if m["isint"].any():
        r = milp(m["c"], constraints=LinearConstraint(A, m["rlo"], m["rup"]), integrality=m["isint"].astype(int),
                 bounds=Bounds(m["lo"], m["up"]), options=dict(mip_rel_gap=1e-9, time_limit=120))
    else:
        eq = m["rlo"] == m["rup"]; le = ~eq & np.isfinite(m["rup"]); ge = ~eq & np.isfinite(m["rlo"])
        Aub = np.vstack([m["A"][le], -m["A"][ge]]); bub = np.concatenate([m["rup"][le], -m["rlo"][ge]])
        r = linprog(m["c"], A_ub=Aub if len(bub) else None, b_ub=bub if len(bub) else None, A_eq=m["A"][eq],
                    b_eq=m["rlo"][eq], bounds=list(zip(m["lo"], m["up"])), method="highs")
    st = {0: "optimal", 2: "infeasible", 3: "unbounded"}.get(r.status, "other")
    return st, (float(r.fun) if st == "optimal" else None)


def violation(m, x):
    ax = m["A"] @ x
    v = max(0.0, float(np.max(m["rlo"] - ax, initial=0)), float(np.max(ax - m["rup"], initial=0)),
            float(np.max(m["lo"] - x, initial=0)), float(np.max(x - m["up"], initial=0)))
    if m["isint"].any():
        v = max(v, float(np.max(np.abs(x[m["isint"]] - np.round(x[m["isint"]])), initial=0)))
    return v


def kpis(m, x):
    val = dict(zip(m["cols"], x))
    s = lambda pre: float(sum(v for k, v in val.items() if k.startswith(pre)))
    return dict(profit=-float(m["c"] @ x), 
                crude_bought_bbl=s("BUY_SWEET") + s("BUY_SOUR"), sweet_bbl=s("BUY_SWEET"), sour_bbl=s("BUY_SOUR"),
                crude_run_bbl=s("RUN_"), gasoline_sold=s("SELL_GAS"), diesel_sold=s("SELL_DSL"),
                fuel_oil_sold=s("SELL_FO"), blendstock_bought=s("BU_G_") + s("HD_D_"),
                unit_on_periods=s("ON_"), startups=s("UP_"))


def run_engine(engine, path):
    tmp = tempfile.mkdtemp(prefix="refstress-")
    js, sol = os.path.join(tmp, "o.json"), os.path.join(tmp, "o.sol")
    for p in (js, sol):
        if os.path.exists(p):
            os.remove(p)
    subprocess.run([engine, path, "--time-limit", "100", "--json", js, "--sol", sol], capture_output=True, timeout=600)
    r = json.load(open(js))
    vals = dict(l.split() for l in open(sol)) if os.path.exists(sol) else None
    return r, vals


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="./taral"); ap.add_argument("--dir", default=os.path.join(here, "refinery_stress"))
    ap.add_argument("--write"); ap.add_argument("--report")
    a = ap.parse_args()
    if a.write:
        os.makedirs(a.write, exist_ok=True)
        for n, sc in scenarios().items():
            write_mps(n, build(sc), os.path.join(a.write, n + ".mps"))
        print("wrote", len(scenarios()), "cases to", a.write)
        return 0
    wrong, report = 0, {}
    print("SYNTHETIC data; linear models only (no nonlinear pooling).")
    for n, sc in scenarios().items():
        path = os.path.join(a.dir, n + ".mps")
        m = parse_mps(path)
        hst, hobj = highs(m)
        truth = sc.get("truth")
        r, vals = run_engine(a.engine, path)
        st = r["status"].lower()
        verdict, why, rep = "ok", "", None
        if truth and hst != truth:
            verdict, why = "wrong", f"oracle says {hst}, construction says {truth}"
        has = vals is not None and (st == "optimal")
        if has:
            if any(c not in vals for c in m["cols"]):
                wrong += 1
                report[n] = dict(verdict="wrong", why="solution file misses columns")
                print(f"{n:28s} WRONG solution file misses columns")
                continue
            x = np.array([float(vals[c]) for c in m["cols"]])
            v = violation(m, x); obj = float(m["c"] @ x)
            if v > 1e-6:
                verdict, why = "wrong", f"point violates model by {v:.3g}"
            elif hst != "optimal":
                verdict, why = "wrong", f"engine optimal but oracle {hst}"
            elif abs(obj - hobj) > 1e-6 * max(1, abs(hobj)):
                verdict, why = "wrong", f"objective {obj:.9g} vs oracle {hobj:.9g}"
            elif abs(obj - r["objective"]) > 1e-6 * max(1, abs(obj)):
                verdict, why = "wrong", "reported objective != objective recomputed from point"
            else:
                rep = kpis(m, x)
        elif hst == "optimal":
            verdict, why = "wrong", f"engine says {st}, oracle optimal {hobj:.9g}"
        elif st not in ("infeasible", "infeasible_or_unbounded") and hst == "infeasible":
            verdict, why = "wrong", f"engine says {st}, oracle infeasible"
        wrong += verdict == "wrong"
        report[n] = dict(verdict=verdict, why=why, engine_status=st, oracle_status=hst, oracle_objective=hobj,
                         engine_objective=r.get("objective"), kpis=rep, synthetic=True)
        k = (" profit %.1f crude_run %.0f" % (rep["profit"], rep["crude_run_bbl"])) if rep else ""
        print(f"{n:28s} {verdict.upper():5s} engine={st:10s} oracle={hst:10s} obj={r.get('objective')}{k} {why}")
    if a.report:
        json.dump(report, open(a.report, "w"), indent=1)
    print("WRONG cases:", wrong)
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main())
