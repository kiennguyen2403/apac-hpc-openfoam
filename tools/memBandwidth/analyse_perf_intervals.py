#!/usr/bin/env python3
"""Memory-stall breakdown of the solver's time-step window from per-rank `perf stat -I 1000 -x,` files.

usage: analyse_perf_intervals.py PERF_DIR SOLVER_LOG STREAM_GBPS

The time-step window is taken from the solver log: from the end of step 1 to the end of the last step
(ClockTime lines, wall clock; perf's interval time starts at process start, slightly before the solver
clock, so margins are dropped at each end). Counts are summed over the window
(perf scales multiplexed counts). Per rank and over all ranks:
  stall share       CYCLE_ACTIVITY.STALLS_TOTAL / cycles
  L2-miss / L3-miss stall share   STALLS_L2_MISS, STALLS_L3_MISS / cycles
  bandwidth-limited OFFCORE_REQUESTS_OUTSTANDING.DATA_RD with >= 4 outstanding (c4) / cycles
  latency-limited   (cycles with >= 1 outstanding data read - c4) / cycles
  page-walk share   DTLB_LOAD_MISSES.WALK_ACTIVE / cycles
  L2 fill traffic   L2_LINES_IN.ALL * 64 B, GB/s summed over ranks, vs the STREAM reference
Caveat: MPI busy-wait inside the window is counted too (cycles that rarely stall on memory), so all
shares are diluted by the MPI fraction of the window (about 10-15% at 1 node).
"""
import re, sys
from pathlib import Path

perf_dir, log, stream = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])

times = [float(m) for m in re.findall(r'ClockTime = ([\d.]+) s', log.read_text())]
if len(times) < 3:
    sys.exit('not enough ClockTime lines')
# perf time = ClockTime + (process start -> Time object, ~1-3 s): skip 2 s after step 1, 1 s before the end
t_lo, t_hi = times[0] + 2.0, times[-1] - 1.0
print(f'time-step window: {t_lo:.1f} s .. {t_hi:.1f} s ({len(times) - 1} steps after step 1)')

NAMES = ['cycles', 'instructions', 'stalls_total', 'stalls_l2', 'stalls_l3', 'rd_any', 'rd_c4',
         'l2_lines_in', 'l3_miss_loads', 'dtlb_walk']


def load(f):
    tot = {k: 0.0 for k in NAMES}
    tmin, tmax = None, None
    for ln in f.read_text().splitlines():
        if not ln or ln.startswith('#'):
            continue
        x = ln.split(',')
        if len(x) < 4:
            continue
        try:
            t = float(x[0])
            v = float(x[1])
        except ValueError:
            continue
        ev = x[3].replace(':u', '')
        if ev not in tot or not (t_lo <= t <= t_hi):
            continue
        tot[ev] += v
        tmin = t if tmin is None else min(tmin, t)
        tmax = t if tmax is None else max(tmax, t)
    tot['span'] = (tmax - tmin + 1.0) if tmin is not None else 0.0
    return tot


ranks = {}
for f in sorted(perf_dir.glob('rank*.txt')):
    ranks[int(f.stem[4:])] = load(f)
R = [r for r in ranks if ranks[r]['cycles'] > 0]
S = {k: sum(ranks[r][k] for r in R) for k in NAMES}
span = max(ranks[r]['span'] for r in R)


def share(k, d=S):
    return 100.0 * d[k] / d['cycles']


print(f'{len(R)} ranks with data, window span {span:.0f} s')
print(f'IPC {S["instructions"] / S["cycles"]:.2f}')
print(f'stalled cycles (any reason)        {share("stalls_total"):5.1f}%')
print(f'  stalled with an L2 miss pending  {share("stalls_l2"):5.1f}%')
print(f'  stalled with an L3 miss pending  {share("stalls_l3"):5.1f}%')
print(f'cycles with >=1 data read to uncore outstanding {share("rd_any"):5.1f}%')
print(f'  bandwidth-limited (>=4 outstanding)         {share("rd_c4"):5.1f}%')
print(f'  latency-limited (1-3 outstanding)           {share("rd_any") - share("rd_c4"):5.1f}%')
print(f'DTLB load page-walk active cycles   {share("dtlb_walk"):5.1f}%')
gbps = S['l2_lines_in'] * 64 / span / 1e9
print(f'L2 fill traffic (from L3 + DRAM), node total: {gbps:.0f} GB/s = {100 * gbps / stream:.0f}% of the STREAM triad '
      f'reference {stream:.0f} GB/s (L3 hits included, so an upper bound for DRAM traffic)')
print(f'demand loads missing L3: {S["l3_miss_loads"] / span / 1e6:.0f} M/s node total '
      f'({S["l3_miss_loads"] * 64 / span / 1e9:.1f} GB/s demand-only, prefetches excluded)')
vals = sorted((share('stalls_l3', ranks[r]), r) for r in R)
print(f'per rank L3-miss stall share: min {vals[0][0]:.1f}% (rank {vals[0][1]}), '
      f'median {vals[len(vals) // 2][0]:.1f}%, max {vals[-1][0]:.1f}% (rank {vals[-1][1]})')
