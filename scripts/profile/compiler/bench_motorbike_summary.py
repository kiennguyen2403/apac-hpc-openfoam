#!/usr/bin/env python3
"""Summarise bench_motorbike_spr.sh: contest timer, perf counters, perf hotspots, residual agreement."""
import csv
import re
import sys
from pathlib import Path

VARIANTS = ('vanilla', 'gcc', 'clang', 'icx')
TIMER = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
RESID = re.compile(r'^(\S+):\s+Solving for (\S+), Initial residual = ([0-9.eE+-]+), Final residual = [0-9.eE+-]+, No Iterations (\d+)')
FP_WIDTH = {'fp_scalar_double': 1, 'fp_128b_packed_double': 2, 'fp_256b_packed_double': 4, 'fp_512b_packed_double': 8}


def timer(log):
    found = TIMER.findall(log.read_text())
    if not found:
        raise SystemExit('missing contest timer in ' + str(log))
    return float(found[-1])


def residuals(log):
    """First-solve initial residual and total linear iterations per field, per time step."""
    steps, current = [], None
    for line in log.read_text().splitlines():
        if line.startswith('Time = '):
            current = {}
            steps.append(current)
            continue
        m = RESID.match(line)
        if m and current is not None:
            key = m.group(2)
            if key not in current:
                current[key] = [float(m.group(3)), 0]
            current[key][1] += int(m.group(4))
    return steps


def perf_stat(directory):
    totals, ranks = {}, 0
    for path in sorted(directory.glob('stat.*.csv')):
        ranks += 1
        with path.open() as stream:
            for row in csv.reader(stream):
                if len(row) < 3 or row[0].startswith('#') or not row[0].strip():
                    continue
                try:
                    value = float(row[0])
                except ValueError:
                    continue   # <not counted> / <not supported>
                name = row[2].split('/')[0].replace(':u', '')
                totals[name] = totals.get(name, 0.0) + value
    return ranks, totals


def hotspots(report, limit=12):
    rows = []
    for line in report.read_text().splitlines():
        m = re.match(r'\s+([0-9.]+)%\s+(\S+)\s+\[\.\]\s+(.*)$', line)
        if m:
            rows.append((float(m.group(1)), m.group(2), m.group(3).strip()))
    return rows[:limit]


def main(root, record_ranks):
    record_ranks = record_ranks.split()
    lines = ['# HPC_motorbike Medium compiler comparison', '',
             '17,444,758 cells, 1 normalsr node, 104 MPI ranks (one per core), 40 simpleFoam iterations.',
             'Metric: contest timer "Average wall-clock time per time step" (from the second step on); lower is better.', '']
    timing = {}
    with (root / 'summary.tsv').open('w') as stream:
        out = csv.writer(stream, delimiter='\t')
        out.writerow(['build', 't1_s', 't2_s', 'mean_s', 'speedup_vs_vanilla', 'perf_trial_s', 'perf_overhead_pct'])
        for v in VARIANTS:
            t = [timer(root / f'run_{v}_{r}' / 'simpleFoam.log') for r in ('t1', 't2')]
            timing[v] = (t, sum(t) / 2, timer(root / f'run_{v}_perf' / 'simpleFoam.log'))
        base = timing['vanilla'][1]
        for v in VARIANTS:
            t, mean, p = timing[v]
            out.writerow([v, t[0], t[1], mean, base / mean, p, 100 * (p / mean - 1)])
    lines += ['| Build | Repeat 1 (s) | Repeat 2 (s) | Mean (s/step) | Speed-up vs vanilla |', '|---|---:|---:|---:|---:|']
    for v in sorted(VARIANTS, key=lambda x: timing[x][1]):
        t, mean, _ = timing[v]
        lines.append(f'| {v} | {t[0]:.4f} | {t[1]:.4f} | {mean:.4f} | {base / mean:.3f}x |')

    lines += ['', '## Hardware counters (perf stat, user space, summed over stat ranks)', '',
              'Whole process including start-up and mesh reading. FP counts are double-precision',
              'arithmetic instructions retired; FMA counts twice. Vector share = fraction of FP elements in 512-bit instructions.', '',
              '| Build | Ranks | Instructions (1e12) | IPC | FP elements (1e12) | 512-bit share | 256-bit share | Scalar share | Cache misses (1e9) | Branch misses (1e9) |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for v in VARIANTS:
        ranks, c = perf_stat(root / f'run_{v}_perf')
        if not ranks or 'instructions' not in c:
            lines.append(f'| {v} | {ranks} | perf stat unavailable | | | | | | | |')
            continue
        elements = {k: c.get(k, 0.0) * w for k, w in FP_WIDTH.items()}
        total = sum(elements.values()) or float('nan')
        lines.append('| %s | %d | %.3f | %.2f | %.3f | %.1f%% | %.1f%% | %.1f%% | %.2f | %.2f |' % (
            v, ranks, c['instructions'] / 1e12, c['instructions'] / c['cycles'], total / 1e12,
            100 * elements['fp_512b_packed_double'] / total, 100 * elements['fp_256b_packed_double'] / total,
            100 * elements['fp_scalar_double'] / total, c.get('cache-misses', 0) / 1e9, c.get('branch-misses', 0) / 1e9))

    lines += ['', '## Hotspots (perf record, cycles:u, rank %s)' % record_ranks[0], '']
    for v in VARIANTS:
        report = root / f'run_{v}_perf' / f'report.{record_ranks[0]}.txt'
        lines += [f'### {v}', '', '| % cycles | Library | Symbol |', '|---:|---|---|']
        rows = hotspots(report) if report.exists() else []
        lines += [f'| {p:.1f} | {dso} | `{sym[:90]}` |' for p, dso, sym in rows] or ['| | perf record unavailable | |']
        lines.append('')

    lines += ['## Numerical agreement', '',
              'Largest relative difference of per-step initial residuals versus vanilla repeat 1, and total linear-solver iterations.', '',
              '| Run | Max rel. residual diff | p iterations | U iterations |', '|---|---:|---:|---:|']
    reference = residuals(root / 'run_vanilla_t1' / 'simpleFoam.log')
    for v in VARIANTS:
        for r in ('t1', 't2', 'perf'):
            steps = residuals(root / f'run_{v}_{r}' / 'simpleFoam.log')
            worst = 0.0
            for a, b in zip(reference, steps):
                for k in a:
                    if k in b and a[k][0] > 0:
                        worst = max(worst, abs(a[k][0] - b[k][0]) / a[k][0])
            p_it = sum(s.get('p', [0, 0])[1] for s in steps)
            u_it = sum(s.get('Ux', [0, 0])[1] for s in steps)
            lines.append(f'| {v} {r} | {worst:.2e} | {p_it} | {u_it} |')
    lines += ['', 'vanilla = module OpenFOAM/2512 (GCC 8.5, OpenMPI 4.1.7); gcc/clang/icx = custom v2512 OptSPR builds',
              'with Intel MPI 2021.17.2, so vanilla-vs-custom also differs in MPI library.', '']
    (root / 'REPORT.md').write_text('\n'.join(lines))
    print('REPORT_SUCCESS ' + str(root / 'REPORT.md'))


if __name__ == '__main__':
    main(Path(sys.argv[1]), sys.argv[2])
