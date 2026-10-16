#!/usr/bin/env python3
"""Fit a per-rank cost model on measured per-step compute, to derive per-cell decomposition weights.

Per-rank per-step compute from the imbalance diagnosis (APS elapsed - MPI, 25 minus 5 steps),
averaged over the normal and permuted mappings (same subdomain on different hardware, so placement
noise averages out). Features per rank from decompScore (same decomposition, case416_w05).

usage: fit_cell_cost.py <diag ROOT> <features .ranks> <REPORT.md>
"""
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

IO_RANKS = {0, 104, 208, 312}
COLS = ['cells', 'procFaces', 'bndFaces', 'neighbours', 'offNode', 'internalFaces', 'wallFaces', 'wallCells']


def aps(p):
    o = {}
    for line in Path(p).read_text().splitlines():
        m = re.match(r'\s*(\d{4})\s+([0-9.]+)\s+([0-9.]+)\s+', line)
        if m:
            o[int(m.group(1))] = float(m.group(2)) - float(m.group(3))
    return o


def main():
    root, feats, report = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    c5, c25, cp = (aps(root / f'aps-t-{n}.txt') for n in ('normal5', 'normal25', 'perm25'))
    rows = {int(r[0]): dict(zip(COLS, map(float, r[1:9])))
            for r in (l.split() for l in Path(feats).read_text().splitlines() if l.strip() and not l.startswith('#'))}
    ranks = [r for r in sorted(rows) if r in c5 and r in c25 and r in cp and r not in IO_RANKS]
    y = np.array([((c25[r] - c5[r]) + (cp[r] - c5[r])) / 40 for r in ranks])
    F = {k: np.array([rows[r][k] for r in ranks]) for k in COLS}
    F['otherBndFaces'] = F['bndFaces'] - F['wallFaces']
    F['totalFaces'] = 2 * F['internalFaces'] + F['procFaces'] + F['bndFaces']

    out = ['# Per-cell cost model from measured per-rank compute', '',
           f'{len(ranks)} ranks; per-step compute (mean of normal and permuted mapping): '
           f'median {np.median(y):.4f} s, max {y.max():.4f} s', '',
           '## Feature summary (per rank)', '', '| feature | min | median | max | corr with compute |', '|---|---:|---:|---:|---:|']
    for k in ['cells', 'internalFaces', 'procFaces', 'wallFaces', 'otherBndFaces', 'wallCells', 'neighbours', 'totalFaces']:
        v = F[k]
        out.append(f'| {k} | {v.min():.0f} | {np.median(v):.0f} | {v.max():.0f} | {np.corrcoef(v, y)[0, 1]:+.2f} |')

    cands = ['cells', 'internalFaces', 'procFaces', 'wallFaces', 'otherBndFaces', 'wallCells', 'neighbours']
    fits = []
    for n in range(1, 5):
        for cols in combinations(cands, n):
            X = np.column_stack([np.ones(len(ranks))] + [F[c] for c in cols])
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            r2 = 1 - ((y - X @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            fits.append((r2, cols, coef))
    fits.sort(key=lambda f: -f[0])
    out += ['', '## Best linear fits (per-step compute ~ a + sum b_i * feature_i)', '',
            '| R2 | features | coefficients (us per unit) |', '|---:|---|---|']
    for r2, cols, coef in fits[:12]:
        out.append(f'| {r2:.3f} | {" + ".join(cols)} | a={coef[0]:.3f} s; '
                   + ', '.join(f'{c}={1e6 * b:.3f}' for c, b in zip(cols, coef[1:])) + ' |')

    # domain-level check: does the best model explain why node domain 0 is cheaper?
    best = fits[0]
    X = np.column_stack([np.ones(len(ranks))] + [F[c] for c in best[1]])
    pred = X @ best[2]
    out += ['', f'## Node-domain means, measured vs best model ({" + ".join(best[1])})', '',
            '| domain | measured | predicted |', '|---|---:|---:|']
    rk = np.array(ranks)
    for d in range(4):
        m = (rk // 104) == d
        out.append(f'| {d} | {y[m].mean():.4f} | {pred[m].mean():.4f} |')
    Path(report).write_text('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
