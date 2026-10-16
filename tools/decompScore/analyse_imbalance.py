#!/usr/bin/env python3
"""Load-imbalance diagnosis from diag_imbalance_4node.sh (APS per-rank reports + I_MPI_DEBUG pin maps).

Per-rank compute = APS elapsed - MPI time. Per-step compute:
    normal = (normal25 - normal5) / 20       (startup cancels)
    perm   = (perm25   - normal5) / 20       (same subdomains, so the same startup work)
If slowness follows the subdomain, normal and perm agree rank by rank; if it follows the hardware,
times agree by physical position (node, cpu) instead.

usage: analyse_imbalance.py <job ROOT> <train.ranks> <REPORT.md>
"""
import re
import sys
from pathlib import Path

import numpy as np

IO_RANKS = {0, 104, 208, 312}
FEATURES = ['cells', 'procFaces', 'bndFaces', 'neighbours', 'offNode']


def aps_ranks(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r'\s*(\d{4})\s+([0-9.]+)\s+([0-9.]+)\s+', line)
        if m:
            out[int(m.group(1))] = float(m.group(2)) - float(m.group(3))
    return out


def pin_map(log):
    """rank -> (host, cpu) from I_MPI_DEBUG>=4 output."""
    out = {}
    for line in Path(log).read_text(errors='replace').splitlines():
        m = re.match(r'\[0\] MPI startup\(\):\s+(\d+)\s+\d+\s+(\S+)\s+\{?(\d+)', line)
        if m:
            out[int(m.group(1))] = (m.group(2).split('.')[0], int(m.group(3)))
    return out


def topo(lscpu):
    t = {}
    for line in Path(lscpu).read_text().splitlines():
        if line.startswith('#'):
            continue
        cpu, core, sock, numa = (int(x) for x in line.split(',')[:4])
        t[cpu] = (sock, numa)
    return t


def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def main():
    root, train, report = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    c5, c25, cp = (aps_ranks(root / f'aps-t-{n}.txt') for n in ('normal5', 'normal25', 'perm25'))
    pm_n, pm_p = pin_map(root / 'run_normal25' / 'simpleFoam.log'), pin_map(root / 'run_perm25' / 'simpleFoam.log')
    tp = topo(root / 'lscpu.txt')
    rows = [l.split() for l in Path(train).read_text().splitlines() if l.strip() and not l.startswith('#')]
    feat = {int(r[0]): dict(zip(FEATURES, map(float, r[1:6]))) for r in rows}

    ranks = sorted(r for r in c5 if r in c25 and r in cp and r in pm_n and r in pm_p and r not in IO_RANKS)
    sn = np.array([(c25[r] - c5[r]) / 20 for r in ranks])
    sp = np.array([(cp[r] - c5[r]) / 20 for r in ranks])
    out = ['# Load-imbalance diagnosis, 4 nodes x 104, Bonus setup, scotch (4 104)/(1 0.5)', '',
           f'{len(ranks)} ranks (collated I/O ranks excluded). Per-step compute = APS (elapsed - MPI), 25 minus 5 steps.', '',
           '| Mapping | median s/step | max s/step | slowest / median | ranks > median + 5% |', '|---|---:|---:|---:|---:|']
    for name, s in (('normal', sn), ('permuted', sp)):
        med = np.median(s)
        out.append(f'| {name} | {med:.4f} | {s.max():.4f} | +{100 * (s.max() / med - 1):.1f}% | {int((s > 1.05 * med).sum())} |')

    # --- does slowness follow the subdomain or the hardware position?
    by_rank = corr(sn, sp)
    pos_n = {pm_n[r]: v for r, v in zip(ranks, sn)}
    pos_p = {pm_p[r]: v for r, v in zip(ranks, sp)}
    common = sorted(set(pos_n) & set(pos_p))
    by_pos = corr([pos_n[k] for k in common], [pos_p[k] for k in common]) if len(common) > 10 else float('nan')
    out += ['', '## Subdomain or hardware?', '',
            f'- correlation of per-rank time, normal vs permuted (same subdomain, different hardware): **{by_rank:+.2f}**',
            f'- correlation of per-position time, normal vs permuted (same node+cpu, different subdomain): **{by_pos:+.2f}**',
            '- high by-rank correlation -> the subdomain decides; high by-position -> the hardware decides.']

    # --- placement groupings (normal and permuted)
    def groups(pm, s, key):
        g = {}
        for r, v in zip(ranks, s):
            host, cpu = pm[r]
            sock, numa = tp.get(cpu, (-1, -1))
            k = {'node': host, 'socket': (host, sock), 'numa': (host, numa)}[key]
            g.setdefault(k, []).append(v)
        return {k: (np.mean(v), len(v)) for k, v in g.items()}
    for key in ('node', 'socket', 'numa'):
        out += ['', f'## Mean per-step compute by {key}', '', f'| {key} | normal | permuted |', '|---|---:|---:|']
        gn, gp = groups(pm_n, sn, key), groups(pm_p, sp, key)
        for k in sorted(set(gn) | set(gp), key=str):
            a = f'{gn[k][0]:.4f} ({gn[k][1]})' if k in gn else '-'
            b = f'{gp[k][0]:.4f} ({gp[k][1]})' if k in gp else '-'
            out.append(f'| {k} | {a} | {b} |')

    # --- subdomain features
    X = {k: np.array([feat[r][k] for r in ranks]) for k in FEATURES}
    savg = (sn + sp) / 2
    out += ['', '## Correlation of per-step compute with subdomain features', '',
            '| feature | normal | permuted | mean of both |', '|---|---:|---:|---:|']
    for k in FEATURES:
        out.append(f'| {k} | {corr(X[k], sn):+.2f} | {corr(X[k], sp):+.2f} | {corr(X[k], savg):+.2f} |')
    A = np.column_stack([np.ones(len(ranks))] + [X[k] for k in FEATURES])
    coef, *_ = np.linalg.lstsq(A, savg, rcond=None)
    r2 = 1 - ((savg - A @ coef) ** 2).sum() / ((savg - savg.mean()) ** 2).sum()
    out += ['', f'Linear fit of the mean of both mappings on all features: R2 = {r2:.2f}']

    # --- slowest ranks
    order = np.argsort(-savg)[:12]
    out += ['', '## Slowest ranks (mean of both mappings)', '',
            '| rank | normal | permuted | node/cpu normal | node/cpu perm | cells | procFaces | bndFaces | neighbours |',
            '|---:|---:|---:|---|---|---:|---:|---:|---:|']
    for i in order:
        r = ranks[i]
        f = feat[r]
        out.append(f'| {r} | {sn[i]:.4f} | {sp[i]:.4f} | {pm_n[r][0]}/{pm_n[r][1]} | {pm_p[r][0]}/{pm_p[r][1]} '
                   f'| {int(f["cells"])} | {int(f["procFaces"])} | {int(f["bndFaces"])} | {int(f["neighbours"])} |')
    med = np.median(savg)
    out += ['', f'Median of both mappings: {med:.4f} s/step']
    Path(report).write_text('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
