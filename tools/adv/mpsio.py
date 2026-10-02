"""MPS writer, independent MPS oracle parser and exact original-model checker for the adversarial suite.

Nothing here shares code with src/mps.cpp. The oracle parser implements the MPS conventions that are
common to the IBM/CPLEX/Gurobi/HiGHS readers and documents the points where the format is ambiguous
(see docs/adversarial.md). Row activities and objectives are recomputed in exact rational arithmetic
from the original file and the solution printed by the engine, so the checker itself adds no rounding.
"""
import math
from fractions import Fraction

INF = float("inf")
BIG = 1e20  # |value| >= BIG on a bound is infinity (CPLEX/Gurobi/HiGHS convention; the format itself is silent)


# ----------------------------------------------------------------------------------------------- writer
def num(v):
    """Format a number for the MPS file. Strings pass through so generators can emit '-0', '1D3', ..."""
    if isinstance(v, str):
        return v
    if isinstance(v, int):
        return str(v)
    if v == INF:
        return "1e+30"
    if v == -INF:
        return "-1e+30"
    return repr(float(v))


class Mps:
    """Builder that can emit every quirk the suite needs; render() gives free-format text."""

    def __init__(self, name="ADV"):
        self.name = name
        self.objsense = None  # None | "MAX" | "MIN"
        self.objsense_style = "section"  # section: 'OBJSENSE\n    MAX'; inline: 'OBJSENSE MAX'
        self.rows = []  # (type, name)
        self.cols = {}  # name -> list of (row, value)
        self.col_order = []
        self.int_cols = set()  # columns emitted inside MARKER blocks
        self.rhs = []  # (row, value)
        self.ranges = []  # (row, value)
        self.bounds = []  # (type, col, value or None)
        self.tail = []  # raw lines appended before ENDATA
        self.endata = True

    def row(self, typ, name):
        self.rows.append((typ, name))
        return name

    def col(self, name, integer=False):
        if name not in self.cols:
            self.cols[name] = []
            self.col_order.append(name)
        if integer:
            self.int_cols.add(name)
        return name

    def add(self, col, row, v):
        self.col(col)
        self.cols[col].append((row, v))

    def render(self):
        out = ["NAME          " + self.name]
        if self.objsense:
            if self.objsense_style == "inline":
                out.append("OBJSENSE " + self.objsense)
            else:
                out += ["OBJSENSE", "    " + self.objsense]
        out.append("ROWS")
        for t, r in self.rows:
            out.append(" %s  %s" % (t, r))
        out.append("COLUMNS")
        in_int = False
        mk = 0
        for c in self.col_order:
            want = c in self.int_cols
            if want != in_int:
                out.append("    MARKER%04d  'MARKER'  '%s'" % (mk, "INTORG" if want else "INTEND"))
                mk += 1
                in_int = want
            ent = self.cols[c]
            for i in range(0, len(ent), 2):
                chunk = ent[i:i + 2]
                out.append("    " + c + "  " + "  ".join("%s  %s" % (r, num(v)) for r, v in chunk))
        if in_int:
            out.append("    MARKER%04d  'MARKER'  'INTEND'" % mk)
        if self.rhs:
            out.append("RHS")
            for r, v in self.rhs:
                out.append("    RHS  %s  %s" % (r, num(v)))
        if self.ranges:
            out.append("RANGES")
            for r, v in self.ranges:
                out.append("    RNG  %s  %s" % (r, num(v)))
        if self.bounds:
            out.append("BOUNDS")
            for t, c, v in self.bounds:
                if v is None:
                    out.append(" %s BND  %s" % (t, c))
                else:
                    out.append(" %s BND  %s  %s" % (t, c, num(v)))
        out += self.tail
        if self.endata:
            out.append("ENDATA")
        return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------------------- oracle parser
class OracleError(Exception):
    pass


class Model:
    """Parsed original model. Bounds/ranges already resolved per the conventions in the module doc."""

    def __init__(self):
        self.name = ""
        self.maximize = False
        self.obj_const = 0.0
        self.row_names = []
        self.row_lo = []
        self.row_up = []
        self.col_names = []
        self.cost = []
        self.col_lo = []
        self.col_up = []
        self.is_int = []
        self.entries = []  # per column: dict row_index -> float (duplicates summed in the order listed)


def _f(s):
    try:
        return float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise OracleError("bad number %r" % s)


def parse_mps(path):
    """Free-format (whitespace) reader. Raises OracleError on anything the spec does not allow."""
    with open(path, "r", newline="") as fh:
        text = fh.read()
    sec = None
    obj = None
    row_type = {}
    row_idx = {}
    m = Model()
    col_idx = {}
    rhs, rng = {}, {}
    bnd = {}  # col -> [lo, up, lo_set, ...]
    in_int = False
    rows_order = []
    seen_end = False
    set_names = {"RHS": None, "RANGES": None, "BOUNDS": None}
    sense_pending = False
    for raw in text.split("\n"):
        raw = raw.rstrip("\r")
        if not raw.strip() or raw[0] == "*":
            continue
        if raw[0] not in " \t":
            t = raw.split()
            sec = t[0]
            if sec == "ENDATA":
                seen_end = True
                break
            if sec == "NAME":
                m.name = t[1] if len(t) > 1 else ""
            elif sec == "OBJSENSE":
                if len(t) > 1:
                    m.maximize = t[1].upper() in ("MAX", "MAXIMIZE")
            elif sec not in ("ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS"):
                raise OracleError("section " + sec)
            continue
        t = raw.split()
        if sec == "OBJSENSE":
            m.maximize = t[0].upper() in ("MAX", "MAXIMIZE")
        elif sec == "ROWS":
            if len(t) != 2:
                raise OracleError("ROWS line")
            typ, name = t[0].upper(), t[1]
            if typ == "N":
                if obj is None:
                    obj = name
                else:
                    row_type[name] = "N"  # extra free rows are ignored
            elif typ in "ELG":
                if name in row_idx or name == obj:
                    raise OracleError("duplicate row " + name)
                row_idx[name] = len(rows_order)
                rows_order.append(name)
                row_type[name] = typ
            else:
                raise OracleError("row type")
        elif sec == "COLUMNS":
            if len(t) >= 3 and t[1] == "'MARKER'":
                in_int = t[2] == "'INTORG'"
                continue
            if len(t) not in (3, 5):
                raise OracleError("COLUMNS line")
            c = t[0]
            if c not in col_idx:
                col_idx[c] = len(m.col_names)
                m.col_names.append(c)
                m.cost.append(0.0)
                m.is_int.append(in_int)
                m.entries.append({})
            j = col_idx[c]
            for k in range(1, len(t), 2):
                r, v = t[k], _f(t[k + 1])
                if r == obj:
                    m.cost[j] += v
                elif r in row_idx:
                    m.entries[j][row_idx[r]] = m.entries[j].get(row_idx[r], 0.0) + v
                elif row_type.get(r) == "N":
                    pass
                else:
                    raise OracleError("unknown row " + r)
        elif sec in ("RHS", "RANGES"):
            d = rhs if sec == "RHS" else rng
            if len(t) % 2 == 1:
                if set_names[sec] is None:
                    set_names[sec] = t[0]
                skip = t[0] != set_names[sec]
                body = t[1:]
            else:
                skip = False
                body = t
            if skip:
                continue
            for k in range(0, len(body) - 1, 2):
                r, v = body[k], _f(body[k + 1])
                if r != obj and r not in row_idx and row_type.get(r) != "N":
                    raise OracleError("unknown row " + r)
                d[r] = v
        elif sec == "BOUNDS":
            typ = t[0].upper()
            novalue = typ in ("FR", "MI", "PL", "BV")
            if novalue:
                body = t[1:]
                if len(body) == 2:
                    sname, c = body
                elif len(body) == 1:
                    sname, c = None, body[0]
                elif len(body) == 3:  # value present anyway (BV 1)
                    sname, c = body[0], body[1]
                else:
                    raise OracleError("BOUNDS line")
                v = 0.0
            else:
                if len(t) == 4:
                    sname, c, v = t[1], t[2], _f(t[3])
                elif len(t) == 3:
                    sname, c, v = None, t[1], _f(t[2])
                else:
                    raise OracleError("BOUNDS line")
            if sname is not None:
                if set_names["BOUNDS"] is None:
                    set_names["BOUNDS"] = sname
                if sname != set_names["BOUNDS"]:
                    continue
            if c not in col_idx:
                raise OracleError("BOUNDS on undeclared column " + c)
            j = col_idx[c]
            lo, up = bnd.get(j, [0.0, INF])
            if typ in ("LO", "LI"):
                lo = v
            elif typ in ("UP", "UI"):
                up = v
                if v < 0 and lo == 0.0:
                    lo = -INF  # CPLEX convention; ambiguous in the format, avoided by generators
            elif typ == "FX":
                lo = up = v
            elif typ == "FR":
                lo, up = -INF, INF
            elif typ == "MI":
                lo = -INF
            elif typ == "PL":
                up = INF
            elif typ == "BV":
                lo, up = 0.0, 1.0
            else:
                raise OracleError("bound type " + typ)
            if typ in ("BV", "LI", "UI"):
                m.is_int[j] = True
            bnd[j] = [lo, up]
    if obj is None:
        raise OracleError("no objective row")
    if not seen_end:
        raise OracleError("missing ENDATA")
    m.obj_const = -rhs.get(obj, 0.0)
    for r in rows_order:
        b = rhs.get(r, 0.0)
        R = rng.get(r)
        typ = row_type[r]
        if typ == "E":
            lo = up = b
            if R is not None:
                if R >= 0:
                    up = b + R
                else:
                    lo = b + R
        elif typ == "L":
            up = b
            lo = b - abs(R) if R is not None else -INF
        else:
            lo = b
            up = b + abs(R) if R is not None else INF
        m.row_names.append(r)
        m.row_lo.append(lo)
        m.row_up.append(up)
    for j in range(len(m.col_names)):
        lo, up = bnd.get(j, [0.0, INF])
        m.col_lo.append(-INF if lo <= -BIG else lo)
        m.col_up.append(INF if up >= BIG else up)
    m.row_lo = [-INF if v <= -BIG else v for v in m.row_lo]
    m.row_up = [INF if v >= BIG else v for v in m.row_up]
    return m


# --------------------------------------------------------------------------------------------- checker
def exact_check(m, x):
    """x: list of floats in column order. Returns dict of worst violations (exact arithmetic).

    viol_row_abs / viol_row_rel : max over rows of the distance outside [lo, up] (rel = divided by 1+|bound hit|)
    viol_bnd_abs / viol_bnd_rel : same for column bounds; int_viol: worst |x - round(x)| over integer columns
    objective: c'x + const recomputed exactly (returned as float)
    """
    n = len(m.col_names)
    xf = [Fraction(v) for v in x]
    act = [Fraction(0)] * len(m.row_names)
    obj = Fraction(m.obj_const)
    for j in range(n):
        xj = xf[j]
        if xj == 0:
            continue
        obj += Fraction(m.cost[j]) * xj
        for i, a in m.entries[j].items():
            act[i] += Fraction(a) * xj
    ra = rr = 0.0
    for i, a in enumerate(act):
        lo, up = m.row_lo[i], m.row_up[i]
        v, b = Fraction(0), 0.0
        if lo != -INF and a < Fraction(lo):
            v, b = Fraction(lo) - a, lo
        elif up != INF and a > Fraction(up):
            v, b = a - Fraction(up), up
        if v:
            fv = float(v)
            ra = max(ra, fv)
            rr = max(rr, fv / (1.0 + abs(b)))
    ba = br = 0.0
    iv = 0.0
    for j in range(n):
        lo, up = m.col_lo[j], m.col_up[j]
        v, b = Fraction(0), 0.0
        if lo != -INF and xf[j] < Fraction(lo):
            v, b = Fraction(lo) - xf[j], lo
        elif up != INF and xf[j] > Fraction(up):
            v, b = xf[j] - Fraction(up), up
        if v:
            fv = float(v)
            ba = max(ba, fv)
            br = max(br, fv / (1.0 + abs(b)))
        if m.is_int[j]:
            iv = max(iv, abs(x[j] - round(x[j])))
    return dict(viol_row_abs=ra, viol_row_rel=rr, viol_bnd_abs=ba, viol_bnd_rel=br, int_viol=iv,
                objective=float(obj))


def read_sol(path, m):
    """Engine .sol file: 'name value' per line, columns of the original model in order."""
    vals = {}
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) == 2:
                vals[p[0]] = float(p[1])
    return [vals.get(c, 0.0) for c in m.col_names]
