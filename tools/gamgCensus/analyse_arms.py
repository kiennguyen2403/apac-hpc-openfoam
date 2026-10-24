#!/usr/bin/env python3
"""Compare gamgCensus arms run on the same mapping (e.g. UCX settings).

usage: analyse_arms.py ROOT arm1 arm2 ...

Per arm: step wall from the solver log, per-rank step compute (wall - MPI),
its spread and its dependence on off-node processor faces (census nbr + host
lines), MPI time, and the rusage deltas over the measured steps (user/system
time, minor faults, voluntary/involuntary context switches, thread count).
"""
import sys, re
from pathlib import Path
import numpy as np

root = Path(sys.argv[1])
arms = sys.argv[2:]
IO = {0, 104, 208, 312}


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
                n = int(x[1])
                d['n'] = n
                d['wall'] = float(x[2]) / n
                d['mpi'] = float(x[3]) / n
                d['step'] = d['wall'] - d['mpi']
            elif x[0] == 'rusage':
                d['ru'] = [float(v) for v in x[1:7]]
                d['threads'] = int(x[8])
            elif x[0] == 'nbr':
                d['nbr'].append((int(x[2]), int(x[3])))
            elif x[0] == 'phase':
                d['ph'][x[1]] = float(x[2]) - float(x[3])
        D[r] = d
    return D


def walltime(arm):
    t = None
    for ln in (root / arm / 'simpleFoam.log').read_text().splitlines():
        if 'Average wall-clock time per time step' in ln:
            t = float(ln.split('=')[1])
    return t


res = {a: load(a) for a in arms}
R = sorted(r for r in res[arms[0]] if r not in IO)
D0 = res[arms[0]]
off = np.array([sum(n for p, n in D0[r]['nbr'] if D0[p]['host'] != D0[r]['host']) for r in R])
print(f'{len(R)} ranks; off-node faces: total {off.sum()}, ranks with any {(off > 0).sum()}, max {off.max()}')

hdr = (f'{"arm":>9} {"s/step":>8} {"comp mean":>9} {"max%":>6} {"sd ms":>6} {"corr off":>8} '
       f'{"us/offF":>7} {"mpi ms":>7} | {"usr s/st":>8} {"sys ms/st":>9} {"minflt/st":>9} '
       f'{"ivcsw/st":>8} {"vcsw/st":>8} {"thr":>4} {"corr sys~off":>12}')
print(hdr)
for a in arms:
    D = res[a]
    c = np.array([D[r]['step'] for r in R])
    m = np.array([D[r]['mpi'] for r in R])
    n = D[R[0]]['n']
    ru = np.array([D[r]['ru'] for r in R]) / n
    thr = [D[r]['threads'] for r in R]
    slope = np.polyfit(off, c, 1)[0]
    print(f'{a:>9} {walltime(a):8.4f} {1e3 * c.mean():9.1f} {100 * (c.max() / np.median(c) - 1):6.1f} '
          f'{1e3 * c.std():6.1f} {np.corrcoef(c, off)[0, 1]:+8.3f} {1e6 * slope:7.2f} {1e3 * m.mean():7.1f} | '
          f'{ru[:, 0].mean():8.3f} {1e3 * ru[:, 1].mean():9.2f} {ru[:, 2].mean():9.1f} '
          f'{ru[:, 5].mean():8.1f} {ru[:, 4].mean():8.1f} {max(set(thr), key=thr.count):>4} '
          f'{np.corrcoef(ru[:, 1], off)[0, 1] if ru[:, 1].std() > 0 else 0:+12.3f}')

base = res[arms[0]]
cb = np.array([base[r]['step'] for r in R])
for a in arms[1:]:
    D = res[a]
    c = np.array([D[r]['step'] for r in R])
    d = c - cb
    hi = off > 3000
    print(f'\n{a} - {arms[0]}: mean compute change {1e3 * d.mean():+.1f} ms; ranks with >3000 off-node '
          f'faces {1e3 * d[hi].mean():+.1f} ms, without off-node faces {1e3 * d[off == 0].mean():+.1f} ms; '
          f'corr(change, off) {np.corrcoef(d, off)[0, 1]:+.3f}')
    for k in ['solve:p', 'lap:U', 'lap:k', 'divDevReff', 'grad:p', 'turbulence', 'solve:U']:
        pa = np.array([base[r]['ph'].get(k, 0) / base[r]['n'] for r in R])
        pb = np.array([D[r]['ph'].get(k, 0) / D[r]['n'] for r in R])
        print(f'   {k:>11}: {1e3 * pa.mean():7.1f} -> {1e3 * pb.mean():7.1f} ms; off-node ranks '
              f'{1e3 * (pb - pa)[hi].mean():+6.2f}, others {1e3 * (pb - pa)[off == 0].mean():+6.2f}')
