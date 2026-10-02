#!/usr/bin/env python3
"""Independent LP KKT verifier (Python standard library only).

Usage: python benchmarks/simplex_certificate_check.py MODEL.mps RESULT.json
Reads the original MPS with its own parser; never trusts reported residuals or
basis flags. Arrays use original first-appearance row/column order. Checks x,
row duals, reduced costs, activities, objective, dual objective, residuals and
complementarity. Min: rc=c-A'y, positive multipliers select lower bounds;
max reverses multiplier signs. Uses the first named RHS/RANGES/BOUNDS set,
including TARAL's convention that a negative UP releases default lower zero.
LP only: integer markers/bounds and quadratic sections are rejected, not relaxed.
Exit 0 means all recomputed KKT measures and exported-field checks <= tolerance.
This verifies an optimal point, not infeasibility or unboundedness certificates.
"""
import argparse
import json
import math
from pathlib import Path

INF = float('inf')
ENDPOINT_BAND = 1e-12  # keep equal to kEndpointBand in src/simplex.cpp


class UnsupportedModel(ValueError):
    """Not an LP; do not retry these files as fixed-column MPS."""


def read_model(path):
    def parse(fixed):
        section = None
        rows, columns, costs, bounds, rhs, ranges, sets = {}, {}, {}, {}, {}, {}, {}
        objective, offset, sense = None, 0.0, 1
        def number(value):
            v = float(value.replace('D', 'E').replace('d', 'e'))
            if not math.isfinite(v):
                raise ValueError('nonfinite MPS coefficient')
            return v
        for raw in Path(path).read_text().splitlines():
            if not raw.strip() or raw.startswith('*'):
                continue
            t = raw.split()
            if not raw[0].isspace():
                section = t[0]
                if section == 'ENDATA':
                    break
                if section not in ('NAME', 'OBJSENSE', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS'):
                    raise UnsupportedModel('unsupported section: ' + section)
                if section != 'OBJSENSE' or len(t) < 2:
                    continue
                t = t[1:]
            elif fixed and section in ('ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS'):
                f = lambda a, b: raw[a:b].strip()
                if section == 'ROWS':
                    t = [f(1, 3), f(4, 12)]
                elif section == 'BOUNDS':
                    t = [f(1, 3), f(4, 12), f(14, 22)] + ([f(24, 36)] if f(24, 36) else [])
                else:
                    t = [f(4, 12), f(14, 22), f(24, 36)] + ([f(39, 47), f(49, 61)] if f(39, 47) else [])
            if section == 'OBJSENSE':
                if t[0] not in ('MIN', 'MINIMIZE', 'MAX', 'MAXIMIZE'):
                    raise ValueError('bad objective sense')
                sense = -1 if t[0].startswith('MAX') else 1
            elif section == 'ROWS':
                if len(t) != 2 or t[0] not in ('N', 'E', 'L', 'G'):
                    raise ValueError('bad row')
                if t[0] == 'N' and objective is None:
                    objective = t[1]
                rows[t[1]] = t[0]
            elif section == 'COLUMNS':
                if "'MARKER'" in t:
                    raise UnsupportedModel('integer markers are not LP certificates')
                if len(t) not in (3, 5):
                    raise ValueError('bad column')
                j = t[0]
                col = columns.setdefault(j, {})
                costs.setdefault(j, 0.0)
                for k in range(1, len(t), 2):
                    r, v = t[k], number(t[k + 1])
                    if r not in rows:
                        raise ValueError('unknown row')
                    if r == objective:
                        costs[j] += v
                    elif rows[r] != 'N':
                        col[r] = col.get(r, 0.0) + v
            elif section in ('RHS', 'RANGES'):
                start = 1 if fixed or len(t) % 2 else 0
                if start and t[0]:
                    if sets.setdefault(section, t[0]) != t[0]:
                        continue
                for k in range(start, len(t), 2):
                    r, v = t[k], number(t[k + 1])
                    if r not in rows:
                        raise ValueError('unknown rhs/range row')
                    if section == 'RHS' and r == objective:
                        offset = -v
                    else:
                        (rhs if section == 'RHS' else ranges)[r] = v
            elif section == 'BOUNDS':
                typ = t[0]
                if typ in ('BV', 'LI', 'UI'):
                    raise UnsupportedModel('integer bounds are not LP certificates')
                no_value = typ in ('FR', 'MI', 'PL')
                if fixed or len(t) == 4 or (len(t) == 3 and no_value):
                    name, j = t[1:3]
                    value = t[3] if len(t) == 4 else None
                else:
                    name, j = '', t[1]
                    value = t[2] if len(t) == 3 else None
                if name and sets.setdefault('BOUNDS', name) != name:
                    continue
                if j not in columns:
                    raise ValueError('unknown bound column')
                v = number(value) if value is not None else 0.0
                lo, up = bounds.get(j, (0.0, INF))
                if typ == 'LO': lo = v
                elif typ == 'UP':
                    up = v
                    if v < 0 and lo == 0: lo = -INF
                elif typ == 'FX': lo = up = v
                elif typ == 'FR': lo, up = -INF, INF
                elif typ == 'MI': lo = -INF
                elif typ == 'PL': up = INF
                else: raise ValueError('unsupported bound: ' + typ)
                bounds[j] = (lo, up)
        row_bounds = {}
        for r, kind in rows.items():
            if kind == 'N': continue
            b, width = rhs.get(r, 0.0), ranges.get(r)
            lo, up = (b, b) if kind == 'E' else ((-INF, b) if kind == 'L' else (b, INF))
            if width is not None:
                if kind == 'L': lo = b - abs(width)
                elif kind == 'G': up = b + abs(width)
                elif width >= 0: up = b + width
                else: lo = b + width
            row_bounds[r] = (lo, up)
        return columns, costs, bounds, row_bounds, offset, sense
    try:
        return parse(False)
    except UnsupportedModel:
        raise
    except (ValueError, IndexError):
        return parse(True)


def verify(path, certificate, tol=1e-8):
    cols, costs, bounds, rows, offset, sense = read_model(path)
    def vector(key, size):
        v = certificate[key]
        if not isinstance(v, list) or len(v) != size or any(
                not isinstance(a, (int, float)) or isinstance(a, bool) or not math.isfinite(a) for a in v):
            raise ValueError('invalid ' + key)
        return v
    x = vector('x', len(cols))
    y = vector('row_dual', len(rows))
    claimed_rc = vector('reduced_cost', len(cols))
    claimed_act = vector('row_activity', len(rows))
    act = {r: 0.0 for r in rows}
    ymap = dict(zip(rows, y))
    rc, pobj = [], offset
    for j, v in zip(cols, x):
        pobj += costs[j] * v
        rc.append(costs[j] - sum(a * ymap[r] for r, a in cols[j].items()))
        for r, a in cols[j].items(): act[r] += a * v
    raw_rows, magnitudes, scaled_rows = [], [], []
    for r, (lo, up) in rows.items():
        raw = max(0., lo-act[r], act[r]-up)
        terms = sum(abs(col.get(r,0.)*v) for col,v in zip(cols.values(),x))
        rhs = lo if act[r] < lo else up if act[r] > up else lo if lo == up else 0.
        if lo <= act[r] <= up and lo != up:
            # Same roundoff-continuous endpoint rule as the engine (kEndpointBand in simplex.cpp).
            dlo = act[r]-lo if math.isfinite(lo) else INF
            dup = up-act[r] if math.isfinite(up) else INF
            if dlo <= dup and dlo <= ENDPOINT_BAND*(terms+abs(lo)): rhs = lo
            elif dup < dlo and dup <= ENDPOINT_BAND*(terms+abs(up)): rhs = up
        magnitude = terms + abs(rhs)
        raw_rows.append(raw); magnitudes.append(magnitude)
        scaled_rows.append(raw/magnitude if magnitude else raw)
    pres = dres = rv = bv = comp = 0.0
    dobj = offset
    cnorm = max((abs(v) for v in costs.values()), default=0.0)
    def bound_check(v, multiplier, lo, up):
        nonlocal pres, dres, dobj, comp
        violation = max(lo - v, v - up, 0.0)
        if violation:
            pres = max(pres, violation / (1 + abs(lo if v < lo else up)))
        signed = sense * multiplier
        if signed:
            bound = lo if signed > 0 else up
            if not math.isfinite(bound):
                dres = max(dres, abs(multiplier) / (1 + cnorm))
            else:
                dobj += multiplier * bound
                comp = max(comp, abs(multiplier * (v - bound)))
        return violation
    for r, multiplier in zip(rows, y):
        rv = max(rv, bound_check(act[r], multiplier, *rows[r]))
    for j, v, multiplier in zip(cols, x, rc):
        bv = max(bv, bound_check(v, multiplier, *bounds.get(j, (0.0, INF))))
    measures = dict(objective=pobj, dual_objective=dobj, primal_res=pres, dual_res=dres,
                    gap=abs(pobj-dobj)/(1+abs(pobj)), complementarity=comp/(1+abs(pobj)),
                    max_row_viol=rv, max_bound_viol=bv,
                    max_row_violation_magnitude_scaled=max(scaled_rows,default=0.))
    mismatches = []
    for key, expected in measures.items():
        reported = certificate.get(key)
        if not isinstance(reported, (int, float)) or isinstance(reported, bool) or not math.isfinite(reported):
            mismatches.append(key)
        elif abs(reported-expected) > tol * (1+abs(expected)):
            mismatches.append(key)
    for key, expected, claimed in [('reduced_cost', rc, claimed_rc), ('row_activity', list(act.values()), claimed_act),
                                     ('row_violation_abs',raw_rows,vector('row_violation_abs',len(rows))),
                                     ('row_term_magnitude',magnitudes,vector('row_term_magnitude',len(rows))),
                                     ('row_violation_magnitude_scaled',scaled_rows,vector('row_violation_magnitude_scaled',len(rows)))]:
        if any(not math.isfinite(a) or abs(a-b) > tol*(1+abs(a)) for a, b in zip(expected, claimed)):
            mismatches.append(key)
    quality = 'unknown' if not all(math.isfinite(v) for v in measures.values()) else (
        'pass' if max(pres,dres,measures['gap'],measures['complementarity']) <= tol else 'fail')
    if certificate.get('certificate_quality') != quality: mismatches.append('certificate_quality')
    ok = certificate.get('status') == 'optimal' and not mismatches and all(
        math.isfinite(v) for v in measures.values()) and max(pres, dres, measures['gap'], measures['complementarity']) <= tol
    return dict(pass_certificate=ok, certificate_quality=quality, recomputed={k: (v if math.isfinite(v) else None)
                                                for k, v in measures.items()},
                mismatched_fields=mismatches)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('model'); ap.add_argument('certificate'); ap.add_argument('--tol', type=float, default=1e-8)
    args = ap.parse_args()
    try:
        if not math.isfinite(args.tol) or args.tol <= 0:
            raise ValueError('tolerance must be finite and positive')
        result = verify(args.model, json.loads(Path(args.certificate).read_text()), args.tol)
    except (ValueError, KeyError, IndexError, TypeError, OSError, OverflowError) as exc:
        result = dict(pass_certificate=False, error=str(exc))
    print(json.dumps(result, allow_nan=False))
    raise SystemExit(0 if result['pass_certificate'] else 1)


if __name__ == '__main__':
    main()
