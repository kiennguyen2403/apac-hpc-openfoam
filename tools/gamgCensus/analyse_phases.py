#!/usr/bin/env python3
"""Per-rank step-phase analysis of a gamgCensus run (phase lines).

usage: analyse_phases.py CENSUS_DIR FEATURES [ranks to show...]

Per rank: step compute = (step wall - step MPI) / steps; per phase: self
compute = (self wall - self MPI) / steps. "other" = step compute minus all
phases (work in simpleFoamFused itself: H, A, flux, fused pEqn steps, relax,
field algebra, boundary updates). Variance share of a phase = cov(phase, step
compute) / var(step compute): the shares add up to 1.
"""
import sys, re
from pathlib import Path
import numpy as np

cdir, ffile = Path(sys.argv[1]), Path(sys.argv[2])
show = [int(r) for r in sys.argv[3:]] or [75, 67, 103, 299, 134, 313]
IO = {0, 104, 208, 312}

data = {}
for f in cdir.glob('rank*.txt'):
    lines = f.read_text().splitlines()
    r = int(re.match(r'# rank (\d+)', lines[0])[1])
    d = {'phases': {}}
    for ln in lines:
        x = ln.split()
        if x[0] == 'steps':
            d['n'], d['wall'], d['mpi'] = int(x[1]), float(x[2]), float(x[3])
        elif x[0] == 'phase':
            d['phases'][x[1]] = (float(x[2]), float(x[3]), int(x[4]))
    data[r] = d

R = sorted(r for r in data if r not in IO)
n = data[R[0]]['n']
print(f'{len(R)} ranks (I/O ranks excluded), {n} measured steps')


def group(name):
    if name.startswith('fo:'):
        return 'functionObjects'
    return name


names = sorted({group(p) for r in R for p in data[r]['phases']})
comp = {k: np.zeros(len(R)) for k in names}
mpi = {k: np.zeros(len(R)) for k in names}
for i, r in enumerate(R):
    for p, (w, m, c) in data[r]['phases'].items():
        comp[group(p)][i] += (w - m) / n
        mpi[group(p)][i] += m / n
step = np.array([(data[r]['wall'] - data[r]['mpi']) / n for r in R])
stepw = np.array([data[r]['wall'] / n for r in R])
stepm = np.array([data[r]['mpi'] / n for r in R])
comp['other'] = step - sum(comp[k] for k in names)
mpi['other'] = stepm - sum(mpi[k] for k in names)
names.append('other')

print(f'step wall {stepw.mean():.4f} s (sd {1e3 * stepw.std():.1f} ms); step compute median '
      f'{np.median(step):.4f}, max +{100 * (step.max() / np.median(step) - 1):.1f}%, '
      f'min {100 * (step.min() / np.median(step) - 1):.1f}%, sd {1e3 * step.std():.1f} ms; '
      f'MPI mean {stepm.mean():.4f}')

vs = np.var(step)
rows = []
for k in names:
    c = comp[k]
    share = np.cov(step, c, bias=True)[0, 1] / vs
    rows.append((share, k))
rows.sort(reverse=True)
print('\nPhase compute per step (ms): mean, sd, max-min, variance share of step compute, '
      'corr with step compute; MPI inside the phase (ms) mean')
print(f'{"phase":>24} {"mean":>8} {"sd":>6} {"range":>7} {"share":>6} {"corr":>6} | {"mpi":>7}')
for share, k in rows:
    c = comp[k] * 1e3
    m = mpi[k] * 1e3
    cr = np.corrcoef(step, c)[0, 1] if c.std() > 0 else 0
    print(f'{k:>24} {c.mean():8.2f} {c.std():6.2f} {c.max() - c.min():7.2f} {share:6.3f} '
          f'{cr:+6.2f} | {m.mean():7.2f}')

feat = {}
for ln in ffile.read_text().splitlines():
    if ln.startswith('#') or not ln.strip():
        continue
    x = ln.split()
    feat[int(x[0])] = [float(v) for v in x[1:]]
fn = ['cells', 'procFaces', 'bndFaces', 'neighbours', 'offNode', 'internalF', 'wallFaces']
F = np.array([feat[r][:7] for r in R])
print('\nCorrelation with subdomain features (phases with variance share >= 0.03 or |share| big):')
print(f'{"":>24}' + ''.join(f'{x:>11}' for x in fn))
print(f'{"stepCompute":>24}' + ''.join(f'{np.corrcoef(step, F[:, j])[0, 1]:+11.3f}' for j in range(7)))
for share, k in rows:
    if abs(share) < 0.03:
        continue
    c = comp[k]
    print(f'{k:>24}' + ''.join(f'{np.corrcoef(c, F[:, j])[0, 1]:+11.3f}' for j in range(7)))

# linear fit of step compute on features
A = np.c_[np.ones(len(R)), F[:, [0, 1, 3, 6]]]
co, *_ = np.linalg.lstsq(A, step, rcond=None)
pred = A @ co
print(f'\nstep compute ~ cells, procFaces, neighbours, wallFaces: R2 {1 - np.var(step - pred) / vs:.3f}; '
      f'us per cell {1e6 * co[1]:.2f}, per proc face {1e6 * co[2]:.2f}, per neighbour {1e6 * co[3]:.1f}, '
      f'per wall face {1e6 * co[4]:.2f}')

print('\nSelected ranks (ms per step): step compute, then the top phases')
top = [k for s, k in rows[:9]]
print(f'{"rank":>6} {"step":>7} ' + ' '.join(f'{k[:12]:>12}' for k in top))
for r in show:
    if r not in R:
        continue
    i = R.index(r)
    print(f'{r:>6} {1e3 * step[i]:7.1f} ' + ' '.join(f'{1e3 * comp[k][i]:12.1f}' for k in top))
