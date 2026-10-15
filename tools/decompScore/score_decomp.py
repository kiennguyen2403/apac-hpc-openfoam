#!/usr/bin/env python3
"""Rank candidate decompositions by the predicted time of their slowest rank.

Model (fitted on per-step cost, one sample per rank):
    compute_r = a + b*cells + c*procFaces + d*neighbours + e*offNodeFaces
where compute_r = per-step (wall - MPI) of rank r from two IPM runs of the same case (25 minus
5 steps, so startup cancels; job 181327812), collated I/O ranks excluded. This feature set
orders the three measured decompositions correctly but overstates sizes ~5.6x
(w100: predicted +6.48%, measured +1.16%), so a calibrated column divides by 5.6.
MPI time is almost entirely waiting for the slowest rank, so a decomposition's score is
max_r compute_r. Only the ordering is trusted; the % column is a ranking index.

usage: score_decomp.py <ipm_diff dir> <train.ranks> <ranksDir> <reference-candidate> <REPORT.md>
"""
import re
import sys
from pathlib import Path

import numpy as np

IO_RANKS = {0, 104, 208, 312}
FEATURES = ['cells', 'procFaces', 'neighbours', 'offNode']
CALIBRATION = 5.6


def read_ranks(path):
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith('#')]
    a = np.array(rows, dtype=float)
    return {'rank': a[:, 0].astype(int), 'cells': a[:, 1], 'procFaces': a[:, 2],
            'bndFaces': a[:, 3], 'neighbours': a[:, 4], 'offNode': a[:, 5]}


def ipm_compute(d):
    xml = next(Path(d).glob('*.ipm.xml')).read_text()
    return {int(r): float(w) - float(m) for r, w, m in
            re.findall(r'mpi_rank="(\d+)".*?<perf wtime="([^"]+)"[^>]*mtime="([^"]+)"', xml, flags=re.S)}


def fit(ipm_diff, train_ranks):
    t5, t25 = ipm_compute(Path(ipm_diff) / 'steps5'), ipm_compute(Path(ipm_diff) / 'steps25')
    t = {r: (t25[r] - t5[r]) / 20.0 for r in t25 if r in t5}
    f = read_ranks(train_ranks)
    keep = [i for i, r in enumerate(f['rank']) if r in t and r not in IO_RANKS]
    y = np.array([t[f['rank'][i]] for i in keep])
    X = np.column_stack([np.ones(len(keep))] + [f[k][keep] for k in FEATURES])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - ((y - X @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return coef, r2, len(keep)


def score(coef, f):
    X = np.column_stack([np.ones(len(f['cells']))] + [f[k] for k in FEATURES])
    pred = X @ coef
    return pred.max(), int(f['rank'][pred.argmax()])


def main():
    ipm_diff, train, ranks_dir, ref, report = sys.argv[1:6]
    coef, r2, n = fit(ipm_diff, train)
    rows = []
    for p in sorted(Path(ranks_dir).glob('*.ranks')):
        f = read_ranks(p)
        s, worst = score(coef, f)
        med = np.median(f['cells'])
        rows.append(dict(name=p.stem, score=s, worst=worst, n=len(f['cells']),
                         cellsImb=100 * (f['cells'].max() / med - 1), maxPF=int(f['procFaces'].max()),
                         maxBF=int(f['bndFaces'].max()), maxNb=int(f['neighbours'].max()),
                         medNb=int(np.median(f['neighbours'])), totPF=int(f['procFaces'].sum() / 2),
                         totOff=int(f['offNode'].sum() / 2)))
    base = next((r['score'] for r in rows if r['name'] == ref), None)
    rows.sort(key=lambda r: r['score'])
    out = ['# Decomposition search: candidates ranked by predicted slowest rank', '',
           f'Per-step model fitted on {n} ranks (IPM job 181327812, 25-5 steps, scotch (4 104)/(1 0.5)), R2 = {r2:.2f}:',
           '`compute = ' + f'{coef[0]:.2f} s + ' + ' + '.join(
               f'{1e6 * c:.3f} us/{k}' for k, c in zip(FEATURES, coef[1:])) + '`', '',
           f'"vs ref" = predicted change of the slowest rank vs `{ref}`; "calibrated" = vs ref / {CALIBRATION} '
           '(the model overstates sizes). Worth a paired run only if calibrated <= -0.7%.', '',
           '| # | Candidate | vs ref | calibrated | slowest rank | cells max/median | max proc faces | max bnd faces '
           '| neighbours max (median) | total proc faces | off-node faces |',
           '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for i, r in enumerate(rows, 1):
        rel = f"{100 * (r['score'] / base - 1):+.2f}%" if base else 'n/a'
        cal = f"{100 * (r['score'] / base - 1) / CALIBRATION:+.2f}%" if base else 'n/a'
        out.append(f"| {i} | {r['name']} | {rel} | {cal} | {r['worst']} | +{r['cellsImb']:.2f}% | {r['maxPF']:,} | {r['maxBF']:,} "
                   f"| {r['maxNb']} ({r['medNb']}) | {r['totPF']:,} | {r['totOff']:,} |")
    Path(report).write_text('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
