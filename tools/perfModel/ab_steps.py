#!/usr/bin/env python3
"""Per-step wall-clock statistics of A/B runs, robust to stalls.

usage: ab_steps.py AB_DIR [AB_DIR...]
  AB_DIR holds <arm>_r<n>/ directories, each with postProcessing.tar.gz (or postProcessing/) containing
  wallClockTimeStatistics/*/wallClockTime.dat (cumulative wall clock per time step, steps 2..N).

Per run: the per-step times (differences of the cumulative column), their mean (= the solver's "Average wall-clock
time per time step"), median, and the stall excess = sum over steps of max(0, step - 1.1*median), i.e. time lost to
transient stalls (steps more than 10% above the run's median). Per arm: means of the runs' means and medians, vs the first arm.
"""
import io, re, statistics as st, sys, tarfile
from pathlib import Path


def step_times(run):
    run = Path(run)
    txt = None
    tgz = run / 'postProcessing.tar.gz'
    if tgz.exists():
        with tarfile.open(tgz) as tf:
            for m in tf:
                if m.name.endswith('wallClockTime.dat'):
                    txt = tf.extractfile(m).read().decode()
    else:
        for p in run.glob('postProcessing/wallClockTimeStatistics/*/wallClockTime.dat'):
            txt = p.read_text()
    if txt is None:
        return None
    cum = [float(l.split()[1]) for l in txt.splitlines() if l.strip() and not l.startswith('#')]
    return cum[:1] + [b - a for a, b in zip(cum, cum[1:])]   # cumulative starts at step 2


for ab in sys.argv[1:]:
    ab = Path(ab)
    arms = {}
    for d in sorted(ab.glob('*_r[0-9]*')):
        m = re.match(r'(.+)_r(\d+)$', d.name)
        t = step_times(d)
        if m and t:
            arms.setdefault(m[1], []).append(t)
    if not arms:
        print(f'{ab}: no per-step data'); continue
    print(f'== {ab.name}')
    print(f'{"run":16s} {"steps":>5s} {"mean":>8s} {"median":>8s} {"stall excess":>12s}')
    summary = {}
    for a, runs in arms.items():
        for i, t in enumerate(runs, 1):
            med = st.median(t)
            stall = sum(max(0.0, x - 1.1 * med) for x in t)
            print(f'{a + "_r" + str(i):16s} {len(t):5d} {st.mean(t):8.4f} {med:8.4f} {stall:12.3f}')
        summary[a] = (st.mean(st.mean(t) for t in runs), st.mean(st.median(t) for t in runs))
    base = next(iter(summary))
    print(f'{"arm":16s} {"mean":>8s} {"vs " + base:>9s} {"median":>8s} {"vs " + base:>9s}')
    for a, (mn, md) in summary.items():
        bm, bd = summary[base]
        print(f'{a:16s} {mn:8.4f} {100 * (mn / bm - 1):+8.2f}% {md:8.4f} {100 * (md / bd - 1):+8.2f}%')
