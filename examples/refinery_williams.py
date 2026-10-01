"""Refinery planning LP from the open literature, solved by the C++ engine.

Source: the refinery optimisation problem in H. P. Williams, "Model Building in Mathematical
Programming" (Wiley). Two crude oils are distilled into naphthas, oils and residuum; naphthas can be
reformed, oils cracked and residuum turned into lube oil; the streams are blended into premium and
regular petrol (octane specifications), jet fuel (vapour-pressure specification), fuel oil (fixed
recipe 10:4:3:1) and lube oil, to maximise daily profit. All coefficients below are the textbook's;
units are barrels per day and pounds per barrel. The published optimal profit is 211,365.13 per day,
which this script checks. This is a literature benchmark, not operational data from any refinery.

The script also re-optimises four what-if scenarios (crude supply, cracker capacity, lube-oil demand,
petrol price) to show feasible re-planning when inputs change.

Usage (repository root): g++ -O3 -std=c++17 -o taral src/*.cpp && python examples/refinery_williams.py --engine ./taral
"""
import argparse, json, os, subprocess, tempfile

YIELD = {  # distillation output per barrel of crude
    "CRUDE1": {"LN": 0.10, "MN": 0.20, "HN": 0.20, "LO": 0.12, "HO": 0.20, "R": 0.13},
    "CRUDE2": {"LN": 0.15, "MN": 0.25, "HN": 0.18, "LO": 0.08, "HO": 0.19, "R": 0.12},
}
REFORM = {"LN": 0.60, "MN": 0.52, "HN": 0.45}  # reformed gasoline per barrel of naphtha
CRACK = {"LO": {"CO": 0.68, "CG": 0.28}, "HO": {"CO": 0.75, "CG": 0.20}}  # cracked oil / gasoline per barrel
OCTANE = {"LN": 90, "MN": 80, "HN": 70, "RG": 115, "CG": 105}
VAPOUR = {"LO": 1.0, "HO": 0.6, "CO": 1.5, "R": 0.05}
FUEL_RECIPE = {"LO": 10 / 18, "CO": 4 / 18, "HO": 3 / 18, "R": 1 / 18}  # share of each barrel of fuel oil
BASE = {"crude1_max": 20000, "crude2_max": 30000, "distill_max": 45000, "reform_max": 10000, "crack_max": 8000,
        "lube_min": 500, "lube_max": 1000, "premium_octane": 94, "regular_octane": 84, "premium_vs_regular": 0.4,
        "jet_vapour": 1.0, "price": {"PMF": 7.00, "RMF": 6.00, "JF": 4.00, "FO": 3.50, "LBO": 1.50}}
PUBLISHED_OPTIMUM = 211365.13


def build(p):
    """Return (rows, columns, objective) with rows as name -> (sense, rhs) and columns as name -> {row: coef}."""
    rows, cols = {}, {}

    def add(col, row, v):
        cols.setdefault(col, {})[row] = cols.get(col, {}).get(row, 0) + v

    for c, y in YIELD.items():  # stream balances: production - uses = 0
        for s, v in y.items():
            add(c, "BAL_" + s, v)
    for s in ("LN", "MN", "HN"):
        add(s + "_RF", "BAL_" + s, -1)
        add(s + "_RF", "BAL_RG", REFORM[s])
        add(s + "_RF", "REFORMER", 1)
        for blend in ("PMF", "RMF"):
            add(f"{s}_{blend}", "BAL_" + s, -1)
    for s in ("LO", "HO"):
        add(s + "_CR", "BAL_" + s, -1)
        add(s + "_CR", "CRACKER", 1)
        for out, v in CRACK[s].items():
            add(s + "_CR", "BAL_" + out, v)
    for s in ("RG", "CG"):
        for blend in ("PMF", "RMF"):
            add(f"{s}_{blend}", "BAL_" + s, -1)
    for s in ("LO", "HO", "CO", "R"):  # jet fuel components
        add(s + "_JF", "BAL_" + s, -1)
    for s, share in FUEL_RECIPE.items():  # fuel oil uses a fixed recipe
        add("FO", "BAL_" + s, -share)
    add("R_LBO", "BAL_R", -1)
    add("R_LBO", "LUBE", 0.5)
    for blend in ("PMF", "RMF"):  # blend definitions and octane specifications
        spec = p["premium_octane"] if blend == "PMF" else p["regular_octane"]
        for s in ("LN", "MN", "HN", "RG", "CG"):
            add(f"{s}_{blend}", "DEF_" + blend, 1)
            add(f"{s}_{blend}", "OCT_" + blend, OCTANE[s] - spec)  # sum (octane - spec) * volume >= 0
    for s in ("LO", "HO", "CO", "R"):
        add(s + "_JF", "VAP_JF", VAPOUR[s] - p["jet_vapour"])  # sum (vp - spec) * volume <= 0
    add("PMF", "PREM_SHARE", 1)
    add("RMF", "PREM_SHARE", -p["premium_vs_regular"])
    add("CRUDE1", "DISTILL", 1)
    add("CRUDE2", "DISTILL", 1)
    for s in ("LN", "MN", "HN", "LO", "HO", "R", "RG", "CO", "CG"):
        rows["BAL_" + s] = ("E", 0)
    for blend in ("PMF", "RMF"):
        add(blend, "DEF_" + blend, -1)
        rows["DEF_" + blend] = ("E", 0)
        rows["OCT_" + blend] = ("G", 0)
    add("JF", "DEF_JF", -1)
    for s in ("LO", "HO", "CO", "R"):
        add(s + "_JF", "DEF_JF", 1)
    rows.update({"DEF_JF": ("E", 0), "VAP_JF": ("L", 0), "PREM_SHARE": ("G", 0), "DISTILL": ("L", p["distill_max"]),
                 "REFORMER": ("L", p["reform_max"]), "CRACKER": ("L", p["crack_max"])})
    add("LBO", "LUBE", -1)
    rows["LUBE"] = ("E", 0)
    bounds = {"CRUDE1": (0, p["crude1_max"]), "CRUDE2": (0, p["crude2_max"]), "LBO": (p["lube_min"], p["lube_max"])}
    return rows, cols, p["price"], bounds


def write_mps(path, p):
    rows, cols, price, bounds = build(p)
    out = ["NAME REFINERY", "OBJSENSE", "    MAX", "ROWS", " N PROFIT"] + [f" {s} {r}" for r, (s, _) in rows.items()]
    out.append("COLUMNS")
    for c, entries in cols.items():
        if c in price:
            out.append(f" {c} PROFIT {price[c]!r}")
        out += [f" {c} {r} {v!r}" for r, v in entries.items()]
    out.append("RHS")
    out += [f" RHS {r} {rhs!r}" for r, (_, rhs) in rows.items() if rhs]
    out.append("BOUNDS")
    for c, (lo, up) in bounds.items():
        if lo:
            out.append(f" LO BND {c} {lo!r}")
        out.append(f" UP BND {c} {up!r}")
    out.append("ENDATA")
    open(path, "w").write("\n".join(out) + "\n")


def solve(engine, path):
    js, sol = path + ".json", path + ".sol"
    subprocess.run([engine, path, "--time-limit", "30", "--json", js, "--sol", sol], capture_output=True, check=False)
    r = json.load(open(js))
    x = dict(line.split() for line in open(sol)) if r["status"] == "optimal" else {}
    return r, {k: float(v) for k, v in x.items()}


def highs(path):
    try:
        import highspy
    except ImportError:
        return None
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.readModel(path)
    h.run()
    return h.getInfo().objective_function_value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="./taral")
    ap.add_argument("--write-mps", help="also save the base model as this MPS file")
    a = ap.parse_args()
    scenarios = {
        "base (published data)": {},
        "crude 2 supply cut to 25,000 bbl/day": {"crude2_max": 25000},
        "cracker capacity cut to 6,000 bbl/day": {"crack_max": 6000},
        "lube oil demand at least 900 bbl/day": {"lube_min": 900},
        "premium petrol price 7.00 -> 7.50": {"price": dict(BASE["price"], PMF=7.50)},
    }
    ok = True
    with tempfile.TemporaryDirectory() as d:
        for i, (name, change) in enumerate(scenarios.items()):
            p = dict(BASE, **change)
            path = os.path.join(d, f"s{i}.mps")
            write_mps(path, p)
            if i == 0 and a.write_mps:
                write_mps(a.write_mps, p)
            r, x = solve(a.engine, path)
            ref = highs(path)
            agree = ref is None or (r["objective"] is not None and abs(r["objective"] - ref) <= 1e-6 * max(1, abs(ref)))
            ok &= r["status"] == "optimal" and agree
            plan = {k: round(x.get(k, 0), 1) for k in ("CRUDE1", "CRUDE2", "PMF", "RMF", "JF", "FO", "LBO")}
            print(f"{name}: status {r['status']}, profit {r['objective']:.2f}"
                  + (f" (HiGHS {ref:.2f})" if ref is not None else "") + f"\n    plan {plan}")
            if i == 0:
                matches = r["objective"] is not None and abs(r["objective"] - PUBLISHED_OPTIMUM) < 0.01
                ok &= matches
                print(f"    published optimum {PUBLISHED_OPTIMUM:.2f}: {'reproduced' if matches else 'NOT reproduced'}")
    print("PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
