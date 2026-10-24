#!/usr/bin/env python3
"""Rank decomposition candidates (decompScore .ranks files) by their predicted slowest rank.

usage: score_kahip.py SCORE_DIR TRAIN_CENSUS TRAIN_RANKS [TRAIN_CENSUS2 TRAIN_RANKS2 ...]

Cost model: per-rank step compute (census "steps" line: wall - MPI, per step) of
measured runs with the current setup (UCX zero-copy off), regressed on the
decompScore features of the same partition: cells, internal faces, wall faces,
processor faces, neighbours. Training pairs: census dir + .ranks file of that
partition (rank cells procFaces bndFaces neighbours offNodeFaces internalFaces
wallFaces wallCells). I/O ranks (0 104 208 312) are excluded from the fit.
Each candidate is then summarised by max/avg cells, neighbours, processor
faces and the predicted max and spread of per-rank compute.
"""
import sys, re
from pathlib import Path
import numpy as np

IO = {0, 104, 208, 312}
FEAT = ['cells', 'internalFaces', 'wallFaces', 'procFaces', 'neighbours']


def read_ranks(f):
    rows = {}
    for ln in Path(f).read_text().splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        x = ln.split()
        r = int(x[0])
        v = [float(t) for t in x[1:]]
        rows[r] = {'cells': v[0], 'procFaces': v[1], 'bndFaces': v[2], 'neighbours': v[3],
                   'offNode': v[4], 'internalFaces': v[5] if len(v) > 5 else np.nan,
                   'wallFaces': v[6] if len(v) > 6 else np.nan}
    return rows


def read_census(d):
    out = {}
    for f in Path(d).glob('rank*.txt'):
        r = int(re.match(r'rank(\d+)', f.stem)[1])
        for ln in f.read_text().splitlines():
            x = ln.split()
            if x and x[0] == 'steps':
                out[r] = (float(x[2]) - float(x[3])) / int(x[1])
    return out


score_dir = Path(sys.argv[1])
pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
X, y, grp = [], [], []
for i, (cdir, rfile) in enumerate(pairs):
    c, rk = read_census(cdir), read_ranks(rfile)
    for r in sorted(c):
        if r in IO or r not in rk:
            continue
        X.append([rk[r][k] for k in FEAT])
        y.append(c[r])
        grp.append(i)
X, y, grp = np.array(X), np.array(y), np.array(grp)
# per-partition offset (run-to-run level differences), shared slopes
A = np.c_[X, np.eye(len(pairs))[grp]]
co, *_ = np.linalg.lstsq(A, y, rcond=None)
pred = A @ co
r2 = 1 - np.var(y - pred) / np.var(y)
slopes = co[:len(FEAT)]
print(f'cost model on {len(y)} ranks from {len(pairs)} run(s): R2 {r2:.3f}, residual sd {1e3 * np.std(y - pred):.1f} ms')
print('  ' + ', '.join(f'{k} {1e6 * s:.2f} us' for k, s in zip(FEAT, slopes)))
for i in range(len(pairs)):
    m = grp == i
    print(f'  run {i}: measured max {1e3 * y[m].max():.1f} ms, predicted max {1e3 * pred[m].max():.1f} ms, '
          f'corr {np.corrcoef(y[m], pred[m])[0, 1]:+.3f}')
off = co[len(FEAT):].mean()

print('\n| candidate | cells max/avg | nbrs avg | nbrs max | proc faces avg | proc faces max | '
      'total proc faces | wall faces max | pred max (ms) | pred sd (ms) | pred max vs scotch_ref |')
print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
res = {}
for f in sorted(score_dir.glob('*.ranks')):
    rk = read_ranks(f)
    R = sorted(r for r in rk if r not in IO)
    F = np.array([[rk[r][k] for k in FEAT] for r in R])
    p = F @ slopes + off
    cells = np.array([rk[r]['cells'] for r in rk])
    nb = np.array([rk[r]['neighbours'] for r in rk])
    pf = np.array([rk[r]['procFaces'] for r in rk])
    wf = np.array([rk[r]['wallFaces'] for r in rk])
    res[f.stem] = (cells.max() / cells.mean(), nb.mean(), nb.max(), pf.mean(), pf.max(), pf.sum() / 2,
                   wf.max(), p.max(), p.std())
ref = res.get('scotch_ref', (None,) * 9)[7]
for k, v in sorted(res.items(), key=lambda kv: kv[1][7]):
    rel = f'{100 * (v[7] / ref - 1):+.2f}%' if ref else '-'
    print(f'| {k} | {v[0]:.4f} | {v[1]:.1f} | {v[2]:.0f} | {v[3]:.0f} | {v[4]:.0f} | {v[5]:.0f} | {v[6]:.0f} | '
          f'{1e3 * v[7]:.1f} | {1e3 * v[8]:.1f} | {rel} |')
