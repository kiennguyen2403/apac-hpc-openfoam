#!/usr/bin/env python3
"""Refit the per-rank model on per-step cost only, then check it against measured decompositions.

Target per rank = (wall - MPI) of the 25-step run minus that of the 5-step run (IPM, same job,
same case), i.e. the compute time of 20 steps with mesh read / startup cancelled.
Validation: predicted change of the slowest rank vs the reference decomposition, compared with
paired 4-node measurements:
    w100 (1 0.01)            +1.16%  (job 180782454)
    w5 (5 1)                 +0.85%  (job 181017811)
    sock211 (4 2 52)/(2 1 1) +0.35%  (job 181017811, the real decomposePar partition)
Prints the fit for several feature sets; a model is usable only if it orders these correctly.

usage: refit_validate.py <ipm_diff dir> <train.ranks> <w100.ranks> <w5.ranks> <sock211 decomposePar log>
"""
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

IO_RANKS = {0, 104, 208, 312}
ALL = ['cells', 'procFaces', 'bndFaces', 'neighbours', 'offNode']


def ipm_compute(d):
    xml = next(Path(d).glob('*.ipm.xml')).read_text()
    return {int(r): float(w) - float(m) for r, w, m in
            re.findall(r'mpi_rank="(\d+)".*?<perf wtime="([^"]+)"[^>]*mtime="([^"]+)"', xml, flags=re.S)}


def read_ranks(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith('#')]
    a = np.array(rows, dtype=float)
    return dict(zip(['rank'] + ALL, a.T))


def ranks_from_log(path):
    log = Path(path).read_text()
    blocks = re.split(r'^Processor (\d+)\n', log, flags=re.M)
    rows = []
    for i in range(1, len(blocks), 2):
        r, b = int(blocks[i]), blocks[i + 1]
        g = lambda k: int(re.search(k + r' = (\d+)', b).group(1))
        off = sum(int(n) for p, n in re.findall(r'shared with processor (\d+) = (\d+)', b) if int(p) // 104 != r // 104)
        rows.append([r, g('Number of cells'), g('Number of processor faces'), g('Number of boundary faces'),
                     g('Number of processor patches'), off])
    return dict(zip(['rank'] + ALL, np.array(rows, dtype=float).T))


def main():
    d, train, w100, w5, sock = sys.argv[1:6]
    t5, t25 = ipm_compute(Path(d) / 'steps5'), ipm_compute(Path(d) / 'steps25')
    f = read_ranks(train)
    keep = [i for i, r in enumerate(f['rank'].astype(int)) if r in t5 and r in t25 and r not in IO_RANKS]
    y = np.array([t25[int(f['rank'][i])] - t5[int(f['rank'][i])] for i in keep]) / 20.0   # s per step
    print(f'per-step compute per rank: min {y.min():.3f} median {np.median(y):.3f} max {y.max():.3f} s '
          f'(slowest/median {100 * (y.max() / np.median(y) - 1):.1f}%)')
    print('corr(feature, per-step compute):', ', '.join(f'{k}={np.corrcoef(f[k][keep], y)[0, 1]:+.2f}' for k in ALL))
    tests = {'w100': (read_ranks(w100), 1.16), 'w5': (read_ranks(w5), 0.85), 'sock211': (ranks_from_log(sock), 0.35)}
    print('\nmeasured vs ref:', ', '.join(f'{k} {v[1]:+.2f}%' for k, v in tests.items()))
    for n in range(1, len(ALL) + 1):
        for cols in combinations(ALL, n):
            X = np.column_stack([np.ones(len(keep))] + [f[c][keep] for c in cols])
            coef, *_ = np.linalg.lstsq(X, y, rcond=None)
            r2 = 1 - ((y - X @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            def score(g): return (np.column_stack([np.ones(len(g['cells']))] + [g[c] for c in cols]) @ coef).max()
            base = score(f)
            pred = {k: 100 * (score(g) / base - 1) for k, (g, _) in tests.items()}
            meas = {k: m for k, (_, m) in tests.items()}
            order_ok = all((pred[a] > pred[b]) == (meas[a] > meas[b]) for a, b in combinations(meas, 2))
            sign_ok = all(p > 0 for p in pred.values())
            print(f"{'+'.join(cols):42s} R2={r2:.2f}  pred " + ' '.join(f'{k} {pred[k]:+.2f}%' for k in tests)
                  + f"  {'ORDER_OK' if order_ok else 'order_wrong'} {'SIGNS_OK' if sign_ok else 'signs_wrong'}")


if __name__ == '__main__':
    main()
