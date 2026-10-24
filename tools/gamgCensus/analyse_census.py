#!/usr/bin/env python3
"""Per-rank / per-level analysis of a gamgCensus run (tools/gamgCensus).

usage: analyse_census.py CENSUS_DIR FEATURES APS_DIR [ranks to show...]
  CENSUS_DIR  directory of rank<N>.txt files
  FEATURES    decompScore features (# rank cells procFaces bndFaces neighbours
              offNodeFaces internalFaces wallFaces wallCells)
  APS_DIR     imbalance diagnosis dir with aps-t-normal25.txt/aps-t-normal5.txt
              (per-step compute = (elapsed - MPI) over 25 minus 5 steps / 20)
"""
import sys, re
from pathlib import Path
import numpy as np

cdir, ffile, apsdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
show = [int(r) for r in sys.argv[4:]] or [75, 103, 299, 134]
IO = {0, 104, 208, 312}
COMP = ['gs', 'amul', 'resid', 'scale', 'ifInit', 'ifUpd', 'rest']
COMM = ['allreduce', 'wait']

ranks = {}
for f in sorted(cdir.glob('rank*.txt')):
    lines = f.read_text().splitlines()
    m = re.match(r'# rank (\d+) gamgSolves (\d+) gamgSolveTime ([\d.]+)', lines[0])
    r, ns, ts = int(m[1]), int(m[2]), float(m[3])
    hdr = lines[1][2:].split()
    lev, rest, out = [], None, None
    for ln in lines[2:]:
        x = ln.split()
        d = dict(zip(hdr, [x[0]] + [float(v) for v in x[1:]]))
        if x[0] == 'level':
            lev.append(d)
        elif x[0] == 'rest':
            rest = d
        else:
            out = d
    lev.sort(key=lambda d: -d['nCells'])
    ranks[r] = dict(nSolves=ns, tSolve=ts, levels=lev, rest=rest, outside=out)

R = sorted(ranks)
ns = ranks[R[0]]['nSolves']
nlev = max(len(ranks[r]['levels']) for r in R)
print(f'{len(R)} ranks, {ns} GAMG solves each, up to {nlev} levels '
      f'(levels 0..{nlev - 1}, last = coarsest, solved by PCG)')


def comp(d):
    return sum(d.get('t_' + k, 0.0) for k in COMP)


def comm(d):
    return sum(d.get('t_' + k, 0.0) for k in COMM)


# per-level table (per solve = per step), ms
print('\nPer level, per step (ms): cells avg/max, max/avg, compute avg/max, max/avg, '
      'comm avg; interface faces avg/max; Allreduce+Wait calls/step')
print(f'{"lev":>3} {"cells avg":>10} {"max":>8} {"ratio":>6} | {"comp avg":>8} '
      f'{"max":>7} {"ratio":>6} | {"comm":>6} | {"ifF avg":>8} {"max":>7} | {"calls":>6}')
lvl_comp = np.zeros((len(R), nlev + 1))
lvl_cells = np.zeros((len(R), nlev))
lvl_if = np.zeros((len(R), nlev))
lvl_comm = np.zeros((len(R), nlev + 1))
for i, r in enumerate(R):
    for l, d in enumerate(ranks[r]['levels']):
        lvl_comp[i, l] = comp(d) / ns
        lvl_comm[i, l] = comm(d) / ns
        lvl_cells[i, l] = d['nCells']
        lvl_if[i, l] = d['nIfaceFaces']
    lvl_comp[i, nlev] = comp(ranks[r]['rest']) / ns
    lvl_comm[i, nlev] = comm(ranks[r]['rest']) / ns
for l in range(nlev + 1):
    c = lvl_comp[:, l] * 1e3
    cm = lvl_comm[:, l] * 1e3
    if l < nlev:
        n = lvl_cells[:, l]
        f = lvl_if[:, l]
        calls = np.mean([(ranks[r]['levels'][l]['n_allreduce'] + ranks[r]['levels'][l]['n_wait'])
                         / ns if l < len(ranks[r]['levels']) else 0 for r in R])
        print(f'{l:>3} {n.mean():10.0f} {n.max():8.0f} {n.max() / n.mean():6.3f} | '
              f'{c.mean():8.2f} {c.max():7.2f} {c.max() / c.mean():6.3f} | {cm.mean():6.2f} | '
              f'{f.mean():8.0f} {f.max():7.0f} | {calls:6.1f}')
    else:
        calls = np.mean([(ranks[r]['rest']['n_allreduce'] + ranks[r]['rest']['n_wait']) / ns
                         for r in R])
        print(f'{"rst":>3} {"":>10} {"":>8} {"":>6} | {c.mean():8.2f} {c.max():7.2f} '
              f'{c.max() / c.mean():6.3f} | {cm.mean():6.2f} | {"":>8} {"":>7} | {calls:6.1f}')

tot_comp = lvl_comp.sum(1)
fine = lvl_comp[:, 0]
coarse = lvl_comp[:, 1:nlev].sum(1)
rest = lvl_comp[:, nlev]
tot_comm = lvl_comm.sum(1)
solve = np.array([ranks[r]['tSolve'] / ns for r in R])
print(f'\nGAMG p solve per step (s): wall avg {solve.mean():.4f}; compute avg {tot_comp.mean():.4f} '
      f'max {tot_comp.max():.4f} (+{100 * (tot_comp.max() / np.median(tot_comp) - 1):.1f}% vs median); '
      f'comm avg {tot_comm.mean():.4f}')
print(f'  fine level compute   avg {fine.mean():.4f} max {fine.max():.4f} '
      f'(+{100 * (fine.max() / np.median(fine) - 1):.1f}%)')
print(f'  coarse levels compute avg {coarse.mean():.4f} max {coarse.max():.4f} '
      f'(+{100 * (coarse.max() / np.median(coarse) - 1):.1f}%)')
print(f'  rest compute          avg {rest.mean():.4f} max {rest.max():.4f}')
out = np.array([comm(ranks[r]['outside']) / ns for r in R])
print(f'  MPI outside GAMG per step: avg {out.mean():.4f}')

# features and APS compute
feat = {}
for ln in ffile.read_text().splitlines():
    if ln.startswith('#') or not ln.strip():
        continue
    x = ln.split()
    feat[int(x[0])] = [float(v) for v in x[1:]]
fnames = ['cells', 'procFaces', 'bndFaces', 'neighbours', 'offNodeFaces', 'internalFaces',
          'wallFaces', 'wallCells']


def aps(fn):
    d = {}
    for ln in (apsdir / fn).read_text().splitlines():
        x = ln.split()
        if len(x) == 6 and x[0].isdigit():
            d[int(x[0])] = float(x[1]) - float(x[2])
    return d


a25, a5 = aps('aps-t-normal25.txt'), aps('aps-t-normal5.txt')
apsc = {r: (a25[r] - a5[r]) / 20 for r in a25 if r in a5}

keep = [i for i, r in enumerate(R) if r not in IO and r in apsc and r in feat]
Rk = [R[i] for i in keep]
y = np.array([apsc[r] for r in Rk])


def corr(a, b):
    return np.corrcoef(a, b)[0, 1]


print(f'\n{len(Rk)} ranks with APS compute (I/O ranks excluded): median {np.median(y):.4f} s/step, '
      f'max +{100 * (y.max() / np.median(y) - 1):.1f}%, sd {100 * y.std() / y.mean():.2f}%')
gc = tot_comp[keep]
print(f'  corr(APS step compute, GAMG compute) = {corr(y, gc):+.3f}; '
      f'sd of GAMG compute {1e3 * gc.std():.1f} ms vs sd of step compute {1e3 * y.std():.1f} ms')
other = y - gc
print(f'  step compute minus GAMG compute: sd {1e3 * other.std():.1f} ms, '
      f'corr with step compute {corr(y, other):+.3f}')
ally = np.var(y)
print(f'  variance share: GAMG {np.cov(y, gc)[0, 1] / ally:.2f}, rest of step {np.cov(y, other)[0, 1] / ally:.2f}')

print('\nCorrelations with subdomain features (ranks above):')
F = np.array([feat[r] for r in Rk])
targets = [('stepCompute', y), ('gamgComp', gc), ('fine', fine[keep]), ('coarse', coarse[keep]),
           ('nonGamg', other),
           ('coarseCells', lvl_cells[keep, 1:].sum(1)), ('coarseIfFaces', lvl_if[keep, 1:].sum(1))]
print(f'{"":>14}' + ''.join(f'{n[:11]:>12}' for n in fnames))
for name, t in targets:
    print(f'{name:>14}' + ''.join(f'{corr(t, F[:, j]):+12.3f}' for j in range(F.shape[1])))
print(f'\ncorr(coarse compute, coarse cells) = {corr(coarse[keep], lvl_cells[keep, 1:].sum(1)):+.3f}; '
      f'corr(coarse compute, coarse interface faces) = {corr(coarse[keep], lvl_if[keep, 1:].sum(1)):+.3f}')
for l in range(1, min(nlev, 6)):
    print(f'  level {l}: corr(compute, cells) {corr(lvl_comp[keep, l], lvl_cells[keep, l]):+.3f}, '
          f'corr(compute, ifFaces) {corr(lvl_comp[keep, l], lvl_if[keep, l]):+.3f}, '
          f'corr(cells, fine cells) {corr(lvl_cells[keep, l], F[:, 0]):+.3f}, '
          f'corr(cells, wallFaces) {corr(lvl_cells[keep, l], F[:, 6]):+.3f}')

print('\nSelected ranks, per step (ms): step compute (APS), GAMG compute fine/coarse/rest, '
      'level cells 0..4, interface faces 0..2')
for r in show:
    i = R.index(r)
    lc = ' '.join(f'{lvl_cells[i, l]:.0f}' for l in range(min(5, nlev)))
    lf = ' '.join(f'{lvl_if[i, l]:.0f}' for l in range(min(3, nlev)))
    print(f'  rank {r:>3}: {1e3 * apsc.get(r, float("nan")):7.1f} | {1e3 * fine[i]:6.1f} '
          f'{1e3 * coarse[i]:6.1f} {1e3 * rest[i]:5.1f} | cells {lc} | ifF {lf}')
