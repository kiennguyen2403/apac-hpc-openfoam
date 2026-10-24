#!/usr/bin/env python3
"""Score candidate decompositions (decompScore .ranks files) with the per-rank cost model fitted on
the imbalance diagnosis of the current setup (job 182137769, FIT.md).

Model 4: compute = -0.284 s + 2.533 us*internalFaces + 5.227 us*procFaces + 4.099 us*wallFaces
                   + 1286.457 us*neighbours                                   (R2 0.472)
Model 3: compute = -0.379 s + 2.724 us*internalFaces + 6.468 us*procFaces + 4.358 us*wallFaces (R2 0.450)
The step is set by the slowest rank, so the score is the predicted maximum (and its excess over the
median). Collated I/O ranks are included (they compute like the others).

usage: score_weights.py <dir with .ranks> [names...]
"""
import sys
from pathlib import Path

MODELS = {
    'm4': (-0.284, {'internalFaces': 2.533e-6, 'procFaces': 5.227e-6, 'wallFaces': 4.099e-6, 'neighbours': 1286.457e-6}),
    'm3': (-0.379, {'internalFaces': 2.724e-6, 'procFaces': 6.468e-6, 'wallFaces': 4.358e-6}),
}
COLS = ['rank', 'cells', 'procFaces', 'bndFaces', 'neighbours', 'offNode', 'internalFaces', 'wallFaces', 'wallCells']


def load(p):
    rows = []
    for line in Path(p).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        rows.append(dict(zip(COLS, (int(x) for x in line.split()))))
    return rows


def main():
    d = Path(sys.argv[1])
    names = sys.argv[2:] or sorted(p.stem for p in d.glob('*.ranks'))
    base = None
    print('| candidate | max cells | max procFaces | total procFaces | max neighbours | '
          'm4 max (s) | m4 max/median | m3 max (s) | m3 max/median | m4 max vs unweighted |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for n in names:
        rows = load(d / f'{n}.ranks')
        res = {}
        for m, (a, coef) in MODELS.items():
            pred = sorted(a + sum(c * r[f] for f, c in coef.items()) for r in rows)
            res[m] = (pred[-1], pred[len(pred) // 2])
        if n == 'unweighted':
            base = res['m4'][0]
        rel = '%+.2f%%' % (100 * (res['m4'][0] / base - 1)) if base else '-'
        print('| %s | %d | %d | %d | %d | %.4f | %+.1f%% | %.4f | %+.1f%% | %s |' % (
            n, max(r['cells'] for r in rows), max(r['procFaces'] for r in rows),
            sum(r['procFaces'] for r in rows) // 2, max(r['neighbours'] for r in rows),
            res['m4'][0], 100 * (res['m4'][0] / res['m4'][1] - 1),
            res['m3'][0], 100 * (res['m3'][0] / res['m3'][1] - 1), rel))


if __name__ == '__main__':
    main()
