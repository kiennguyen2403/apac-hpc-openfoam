#!/usr/bin/env python3
"""Store-side memory traffic of the solver's time-step window from per-rank `perf stat -I 1000 -x,` files.

usage: analyse_rfo.py PERF_DIR SOLVER_LOG

Events (user mode, summed over the window steps 2..N and over ranks):
  l2_lines_in   L2_LINES_IN.ALL          every line filled into L2 (demand loads, RFOs, prefetches)
  rfo_all       L2_RQSTS.ALL_RFO          store (read-for-ownership) requests reaching L2
  rfo_miss      L2_RQSTS.RFO_MISS         ... of those, missing L2 (line fetched from L3/DRAM)
  offcore_rfo   OFFCORE_REQUESTS.DEMAND_RFO  demand RFOs sent beyond L2
  hwpf_miss     L2_RQSTS.HWPF_MISS        L2 hardware-prefetch requests missing L2 (includes prefetches
                                          issued ahead of store streams)
  wb_dirty      L2_LINES_OUT.NON_SILENT   modified lines written back from L2
What streaming stores can save is the fill of lines that are only written. Upper bound: dirty
write-backs (every written line was filled once; lines that are read and then written need that fill
anyway). Demand RFO misses are a lower bound (store streams the prefetcher covers appear as prefetch fills).
"""
import re, sys
from pathlib import Path

pdir, log = Path(sys.argv[1]), Path(sys.argv[2])
clock = [float(x) for x in re.findall(r'ClockTime = ([\d.]+) s', log.read_text())]
t_lo, t_hi = clock[0] + 2.0, clock[-1] - 1.0
EV = ['cycles', 'instructions', 'l2_lines_in', 'rfo_all', 'rfo_miss', 'offcore_rfo', 'hwpf_miss', 'wb_dirty']
S = {e: 0.0 for e in EV}
span = 0.0
nr = 0
for f in pdir.glob('rank*.txt'):
    ts = []
    for ln in f.read_text().splitlines():
        x = ln.split(',')
        if len(x) < 4 or ln.startswith('#'):
            continue
        try:
            t, v = float(x[0]), float(x[1])
        except ValueError:
            continue
        e = x[3].replace(':u', '')
        if e in S and t_lo <= t <= t_hi:
            S[e] += v
            ts.append(t)
    if ts:
        nr += 1
        span = max(span, max(ts) - min(ts) + 1.0)
steps = len(clock) - 1
step_s = (clock[-1] - clock[0]) / steps
print(f'{nr} ranks, window {span:.0f} s, {step_s:.2f} s/step')
for e in EV:
    if S[e] == 0:
        print(f'WARNING: {e} counted zero (encoding not valid on this CPU?)')
fill = S['l2_lines_in']
gb = lambda n: n * 64 / span / 1e9
print(f'L2 fills (all)                 {gb(fill):7.1f} GB/s  (100%)')
print(f'demand RFO misses (L2)         {gb(S["rfo_miss"]):7.1f} GB/s  ({100 * S["rfo_miss"] / fill:5.1f}% of fills)  <- lower bound')
print(f'demand RFOs sent beyond L2     {gb(S["offcore_rfo"]):7.1f} GB/s  ({100 * S["offcore_rfo"] / fill:5.1f}% of fills)')
print(f'L2 prefetch misses             {gb(S["hwpf_miss"]):7.1f} GB/s  ({100 * S["hwpf_miss"] / fill:5.1f}% of fills)')
print(f'dirty write-backs              {gb(S["wb_dirty"]):7.1f} GB/s  ({100 * S["wb_dirty"] / fill:5.1f}% of fills)  <- upper bound')
if S['rfo_all']:
    print(f'RFO requests hitting L2        {100 * (1 - S["rfo_miss"] / S["rfo_all"]):5.1f}%')
print(f'IPC {S["instructions"] / S["cycles"]:.2f}')
