#!/usr/bin/env python3
"""Compare per-rank compute between two rank-to-node mappings (gamgCensus runs).

usage: analyse_mapping.py ROOT [armA armB]   (default arms: normal rot52)

Per arm and rank: step compute = (step wall - step MPI) / steps, host (census
"host" line), off-node processor faces (census "nbr" lines + that arm's
rank->host map). Tests whether the per-rank compute change between arms
follows the change in off-node faces.
"""
import sys, re
from pathlib import Path
import numpy as np

root = Path(sys.argv[1])
arms = sys.argv[2:4] or ['normal', 'rot52']


def load(arm):
    D = {}
    for f in (root / arm / 'census').glob('rank*.txt'):
        r = int(re.match(r'rank(\d+)', f.stem)[1])
        d = {'nbr': [], 'ph': {}}
        for ln in f.read_text().splitlines():
            x = ln.split()
            if not x:
                continue
            if x[0] == 'host':
                d['host'] = x[1]
            elif x[0] == 'steps':
                d['n'] = int(x[1])
                d['step'] = (float(x[2]) - float(x[3])) / int(x[1])
                d['mpi'] = float(x[3]) / int(x[1])
            elif x[0] == 'nbr':
                d['nbr'].append((int(x[2]), int(x[3])))
            elif x[0] == 'phase':
                d['ph'][x[1]] = (float(x[2]) - float(x[3]))
        D[r] = d
    for d in D.values():
        for k in d['ph']:
            d['ph'][k] /= d['n']
    return D


A, B = load(arms[0]), load(arms[1])
R = sorted(A)


def offnode(D, r):
    return sum(n for p, n in D[r]['nbr'] if D[p]['host'] != D[r]['host'])


def ioranks(D):
    # one collated I/O rank per host: lowest rank of each host block in the
    # normal mapping is 0/104/208/312; exclude those in both arms
    return {0, 104, 208, 312}


io = ioranks(A)
K = [r for r in R if r not in io]
ca = np.array([A[r]['step'] for r in K])
cb = np.array([B[r]['step'] for r in K])
oa = np.array([offnode(A, r) for r in K])
ob = np.array([offnode(B, r) for r in K])
print(f'{len(K)} ranks; mean step compute {arms[0]} {1e3 * ca.mean():.1f} ms, '
      f'{arms[1]} {1e3 * cb.mean():.1f} ms; mean MPI {1e3 * np.mean([A[r]["mpi"] for r in K]):.1f} / '
      f'{1e3 * np.mean([B[r]["mpi"] for r in K]):.1f} ms')
print(f'max step compute: {arms[0]} {1e3 * ca.max():.1f} (+{100 * (ca.max() / np.median(ca) - 1):.1f}%), '
      f'{arms[1]} {1e3 * cb.max():.1f} (+{100 * (cb.max() / np.median(cb) - 1):.1f}%)')
print(f'total off-node faces: {arms[0]} {oa.sum()}, {arms[1]} {ob.sum()}; ranks with any: '
      f'{(oa > 0).sum()} / {(ob > 0).sum()}')
print(f'corr per rank compute {arms[0]} vs {arms[1]}: {np.corrcoef(ca, cb)[0, 1]:+.3f}')
print(f'corr(compute, off-node faces) within arm: {arms[0]} {np.corrcoef(ca, oa)[0, 1]:+.3f}, '
      f'{arms[1]} {np.corrcoef(cb, ob)[0, 1]:+.3f}')
dc, do = cb - ca, ob - oa
print(f'corr(change in compute, change in off-node faces) = {np.corrcoef(dc, do)[0, 1]:+.3f}')
co = np.polyfit(do, dc, 1)
print(f'  slope {1e6 * co[0]:.2f} us per off-node face (intercept {1e3 * co[1]:+.2f} ms)')
gained = do > 2000
lost = do < -2000
print(f'  ranks gaining >2000 off-node faces: {gained.sum()}, mean compute change '
      f'{1e3 * dc[gained].mean():+.1f} ms; losing >2000: {lost.sum()}, '
      f'{1e3 * dc[lost].mean():+.1f} ms; others {1e3 * dc[~gained & ~lost].mean():+.1f} ms')

# which phases move with off-node faces
keys = sorted({k for r in K for k in A[r]['ph']})
rows = []
for k in keys:
    pa = np.array([A[r]['ph'].get(k, 0) for r in K])
    pb = np.array([B[r]['ph'].get(k, 0) for r in K])
    if pa.mean() < 0.002:
        continue
    d = pb - pa
    rows.append((np.corrcoef(d, do)[0, 1] if d.std() > 0 else 0, k, 1e3 * pa.mean(),
                 1e3 * d[gained].mean() if gained.any() else 0,
                 1e3 * d[lost].mean() if lost.any() else 0))
rows.sort(reverse=True)
print('\nphase compute change vs off-node change: corr, phase, mean ms, '
      'change for gaining ranks, for losing ranks (ms)')
for c, k, m, g, l in rows:
    print(f'  {c:+.3f} {k:>20} {m:8.1f} {g:+7.2f} {l:+7.2f}')
