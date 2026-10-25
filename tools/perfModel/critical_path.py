#!/usr/bin/env python3
"""Where a time step goes: average compute, imbalance, communication (critical-path view).

usage: critical_path.py [PHASE_CENSUS_TAR]   (default: diag/phase_census_4node_182140240.tar.gz)

Segments of a step: each census phase (self compute = self wall - self MPI), with the pressure
solve split into its GAMG levels (Gauss-Seidel, Amul, scaling, interface packing; per level) and
the rest of the GAMG solve, plus "other" (step compute not inside any phase). Every segment ends
in a synchronisation (global sums, halo exchanges), so to first order
  step = sum over segments of max over ranks (segment compute) + communication
and it splits into
  average compute                       mean over ranks of the step compute
  whole-step imbalance                  max - mean of the step compute
  per-segment imbalance                 sum of segment maxima - max of the step compute
                                        (different ranks are slowest in different segments)
  communication and finer imbalance     the rest: latency/transfer, waits inside a segment
The census totals are sums over 9 steps per segment, so per-segment maxima of totals understate
cycle-by-cycle maxima: the per-segment imbalance is a lower bound and the last term an upper bound.
"""
import re, sys, tarfile
from statistics import mean

TAR = sys.argv[1] if len(sys.argv) > 1 else \
    '/scratch/ya41/kn3896/openfoam/contest-results/diag/phase_census_4node_182140240.tar.gz'
IO = {0, 104, 208, 312}
GCOMP = ['t_gs', 't_amul', 't_resid', 't_scale', 't_ifInit', 't_ifUpd']
MPI = ['t_allreduce', 't_wait', 't_post', 't_coll']

ranks = {}
with tarfile.open(TAR) as tf:
    for m in tf:
        mm = re.search(r'census/rank(\d+)\.txt$', m.name)
        if not mm:
            continue
        d = {'levels': [], 'phases': {}}
        cols = None
        for ln in tf.extractfile(m).read().decode().splitlines():
            x = ln.split()
            if ln.startswith('# level'):
                cols = x[2:]
            elif x[0] == 'level':
                d['levels'].append(dict(zip(cols, map(float, x[1:]))))
            elif x[0] in ('rest', 'outside'):
                d[x[0]] = dict(zip(cols, map(float, x[1:])))
            elif x[0] == 'phase':
                d['phases'][x[1]] = (float(x[2]), float(x[3]))
            elif x[0] == 'steps':
                d['n'], d['wall'], d['mpi'] = int(x[1]), float(x[2]), float(x[3])
        ranks[int(mm[1])] = d

R = sorted(r for r in ranks if r not in IO)
n = ranks[R[0]]['n']


def segments(d):
    s = {}
    for name, (w, mp) in d['phases'].items():
        s[name] = (w - mp) / n
    g = 0.0
    for i, L in enumerate(d['levels']):
        c = sum(L[k] for k in GCOMP) / n
        s[f'gamg:L{i}'] = c
        g += c
    s['gamg:rest'] = d['rest']['t_rest'] / n
    s['solve:p'] -= g + s['gamg:rest']          # what is left of the p solve outside the levels
    step = (d['wall'] - d['mpi']) / n
    s['other'] = step - sum(s.values())
    return s, step


seg = {r: segments(ranks[r])[0] for r in R}
comp = {r: segments(ranks[r])[1] for r in R}
T = mean(ranks[r]['wall'] for r in R) / n
c_mean = mean(comp.values())
c_max = max(comp.values())
names = sorted(set().union(*[s.keys() for s in seg.values()]))
seg_max = {k: max(seg[r].get(k, 0.0) for r in R) for k in names}
seg_mean = {k: mean(seg[r].get(k, 0.0) for r in R) for k in names}
S = sum(seg_max.values())

ms = lambda x: f'{1000 * x:7.1f} ms  {100 * x / T:5.1f}%'
print(f'{len(R)} ranks (I/O ranks excluded), {n} steps, step {1000 * T:.1f} ms with the census on')
print(f'  average compute                    {ms(c_mean)}')
print(f'  whole-step imbalance (max - mean)  {ms(c_max - c_mean)}')
print(f'  per-segment imbalance (>=)         {ms(S - c_max)}')
print(f'  communication + finer imbalance    {ms(T - S)}')

# MPI time by category, all ranks vs the 20 slowest (by compute)
slow = sorted(R, key=lambda r: -comp[r])[:20]


def mpi_cat(d):
    out = {}
    for k in MPI:
        out[k] = (sum(L[k] for L in d['levels']) + d['rest'][k] + d['outside'][k]) / n
    return out


cat_all = {k: mean(mpi_cat(ranks[r])[k] for r in R) for k in MPI}
cat_slow = {k: mean(mpi_cat(ranks[r])[k] for r in slow) for k in MPI}
print('\nMPI time per step by category      all ranks    20 slowest ranks')
for k in MPI:
    print(f'  {k[2:]:12s}                    {1000 * cat_all[k]:7.1f} ms   {1000 * cat_slow[k]:7.1f} ms')
print(f'  total                           {1000 * sum(cat_all.values()):7.1f} ms   {1000 * sum(cat_slow.values()):7.1f} ms')

print('\nlargest per-segment imbalances (max - mean over ranks, ms/step):')
for k in sorted(names, key=lambda k: -(seg_max[k] - seg_mean[k]))[:12]:
    print(f'  {k:22s} mean {1000 * seg_mean[k]:7.1f}  max {1000 * seg_max[k]:7.1f}  imbalance {1000 * (seg_max[k] - seg_mean[k]):6.1f}')

print('\nceilings (upper bounds on what a perfect fix could save):')
print(f'  perfect balance in every segment   {ms(S - c_mean)}  (decomposition / partitioning work)')
print(f'  zero-cost communication            {ms(T - S)}  (MPI / network tuning; includes finer imbalance)')
