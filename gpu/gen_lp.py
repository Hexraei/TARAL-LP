#!/usr/bin/env python3
"""SYNTHETIC LP generator for the TARAL-LP GPU crossover study (fixed seed; not real-world data).

Families (both feasible and bounded by construction):
  transport  S sources x D sinks, x_ij >= 0, cost = Euclidean distance between random points in the unit square.
             sum_j x_ij <= supply_i,  sum_i x_ij >= demand_j,  total supply = 1.25 x total demand.
             rows S + D, cols S*D, nnz 2*S*D (long rows, 2 entries per column).
  packing    max p'x  s.t.  A x <= b,  0 <= x <= 1,  each column has K entries in random rows, values U(0.1, 1),
             b = A x0 * U(1, 1.5) for a random x0 in [0, 1] (x0 is feasible).  Written as min (-p)'x.
             rows n/2, cols n, nnz K*n.

Usage: python gpu/gen_lp.py transport S D OUT.mps [--seed 1]
       python gpu/gen_lp.py packing N OUT.mps [--k 5] [--seed 1]
Writes free-format MPS.
"""
import argparse
import numpy as np


def fmt(v):
    return np.char.mod('%.12g', v)


def write_mps(path, name, row_types, row_names, cols, rhs, upper=None):
    """cols: list of (col_names array, row_name array, value array) blocks, already grouped by column."""
    with open(path, 'w') as f:
        f.write(f'NAME {name}\nROWS\n N OBJ\n')
        f.write(''.join(f' {t} {r}\n' for t, r in zip(row_types, row_names)))
        f.write('COLUMNS\n')
        for cn, rn, val in cols:
            lines = np.char.add(np.char.add(np.char.add(np.char.add('    ', cn), ' '), np.char.add(rn, ' ')), fmt(val))
            f.write('\n'.join(lines.tolist()))
            f.write('\n')
        f.write('RHS\n')
        rr, rv = rhs
        f.write(''.join(f'    RHS {r} {v:.12g}\n' for r, v in zip(rr, rv)))
        if upper is not None:
            f.write('BOUNDS\n')
            cn, uv = upper
            lines = np.char.add(np.char.add(np.char.add(' UP BND ', cn), ' '), fmt(uv))
            f.write('\n'.join(lines.tolist()))
            f.write('\n')
        f.write('ENDATA\n')


def transport(S, D, path, seed):
    rng = np.random.default_rng(seed)
    ps, pd = rng.random((S, 2)), rng.random((D, 2))
    demand = rng.uniform(1, 10, D)
    supply = rng.uniform(1, 10, S)
    supply *= 1.25 * demand.sum() / supply.sum()
    snames = np.array([f'S{i}' for i in range(S)])
    dnames = np.array([f'D{j}' for j in range(D)])
    blocks = []
    chunk = max(1, 2_000_000 // D)
    for i0 in range(0, S, chunk):
        i1 = min(S, i0 + chunk)
        ii = np.repeat(np.arange(i0, i1), D)
        jj = np.tile(np.arange(D), i1 - i0)
        cost = np.sqrt(((ps[ii] - pd[jj]) ** 2).sum(axis=1))
        cn = np.char.add('X', (ii * D + jj).astype(str))
        # three lines per column: objective, supply row, demand row
        cn3 = np.repeat(cn, 3)
        rn3 = np.stack([np.full(len(ii), 'OBJ'), snames[ii], dnames[jj]], axis=1).ravel()
        v3 = np.stack([cost, np.ones(len(ii)), np.ones(len(ii))], axis=1).ravel()
        blocks.append((cn3, rn3, v3))
    write_mps(path, f'TRANSPORT_{S}x{D}', ['L'] * S + ['G'] * D, list(snames) + list(dnames), blocks,
              (list(snames) + list(dnames), list(supply) + list(demand)))


def packing(N, K, path, seed):
    rng = np.random.default_rng(seed)
    M = max(1, N // 2)
    # K random rows per column (sorted); repeated hits within a column are dropped
    rows = np.sort(rng.integers(0, M, (N, K)), axis=1)
    dup = np.zeros_like(rows, dtype=bool)
    dup[:, 1:] = rows[:, 1:] == rows[:, :-1]
    vals = rng.uniform(0.1, 1.0, (N, K))
    vals[dup] = 0  # duplicate row hits are dropped below
    x0 = rng.random(N)
    ax0 = np.zeros(M)
    np.add.at(ax0, rows.ravel(), (vals * x0[:, None]).ravel())
    b = ax0 * rng.uniform(1.0, 1.5, M) + 0.1
    p = rng.uniform(0.5, 1.5, N)
    rnames = np.array([f'R{i}' for i in range(M)])
    blocks = []
    chunk = 400_000
    for j0 in range(0, N, chunk):
        j1 = min(N, j0 + chunk)
        jj = np.arange(j0, j1)
        keep = ~dup[j0:j1]
        cn = np.char.add('X', jj.astype(str))
        nk = keep.sum(axis=1)
        cn_all = np.concatenate([cn[:, None], np.repeat(cn, K).reshape(-1, K)], axis=1)
        rn_all = np.concatenate([np.full((j1 - j0, 1), 'OBJ'), rnames[rows[j0:j1]]], axis=1)
        v_all = np.concatenate([-p[j0:j1, None], vals[j0:j1]], axis=1)
        mask = np.concatenate([np.ones((j1 - j0, 1), bool), keep], axis=1)
        blocks.append((cn_all[mask], rn_all[mask], v_all[mask]))
    cnames = np.char.add('X', np.arange(N).astype(str))
    write_mps(path, f'PACKING_{N}', ['L'] * M, list(rnames), blocks, (list(rnames), list(b)),
              upper=(cnames, np.ones(N)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family', choices=['transport', 'packing'])
    ap.add_argument('args', nargs='+')
    ap.add_argument('--k', type=int, default=5)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    if a.family == 'transport':
        transport(int(a.args[0]), int(a.args[1]), a.args[2], a.seed)
    else:
        packing(int(a.args[0]), a.k, a.args[1], a.seed)


if __name__ == '__main__':
    main()
