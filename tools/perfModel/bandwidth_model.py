#!/usr/bin/env python3
"""Memory-traffic model of the linear-solver kernels, checked against the measured A/B results.

usage: bandwidth_model.py [PHASE_CENSUS_TAR]   (default: diag/phase_census_4node_182140240.tar.gz)

Model: the solver kernels are memory-bandwidth bound (job 182171696: 72% of the STREAM triad, 72%
of cycles with >= 4 reads outstanding), so a kernel's time is proportional to the bytes it moves.
An optimisation that changes a kernel's bytes from B to B' changes its time by t * (B' - B) / B,
with t the kernel's measured time. Nothing is fitted: the kernel times come from the census
(4 nodes x 104 ranks, Bonus build just before compactAddr, self times with MPI and nested calls
excluded), the bytes from counting the loops, and the predictions are compared with A/B results
that the model never sees.

Bytes per cell (f = faces per cell of the level; 8-byte scalars, 4-byte labels; a write to a line
not yet in cache costs a read for ownership (RFO) plus a write-back):
  Gauss-Seidel call of s sweeps, folded copy (stepReuse):  24 + s(44 + 20f) + 8(s-1)
      copy bPrime = source once (8 + 8 + RFO 8); per sweep ownerStart 4, diag 8, bPrime 8+8,
      psi 8+8; per face upper 8, lower 8, upperAddr 4; later sweeps re-read source 8
  same, stock (copy every sweep):                          s(24 + 44 + 20f)
  Amul:                         32 + r + 24f   diag, psi, Apsi (+RFO); per face lower, upper,
                                               lowerAddr, upperAddr; r = psi/Apsi re-read in the
                                               face loop if evicted, 0..24 (the model's range)
  residual (+ its norm):        48 + r + 24f
  fused Ux/Uy/Uz Gauss-Seidel, 15 sweeps in one call (v3): 72 + 15(108 + 20f) + 24*14
  fused 3-component residual (+ norms):                   128 + 3r + 24f
Compact addressing (uint16 offsets, uint8 faces per cell): -s(3 + 2f) per Gauss-Seidel call,
-(6f - 1) per Amul.
"""
import re, sys, tarfile
from statistics import mean

TAR = sys.argv[1] if len(sys.argv) > 1 else \
    '/scratch/ya41/kn3896/openfoam/contest-results/diag/phase_census_4node_182140240.tar.gz'
IO = {0, 104, 208, 312}
R_LO, R_HI = 0.0, 24.0

# ---------------------------------------------------------------- census
ranks = {}
with tarfile.open(TAR) as tf:
    for m in tf:
        mm = re.search(r'census/rank(\d+)\.txt$', m.name)
        if not mm:
            continue
        r = int(mm[1])
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
                d['phases'][x[1]] = (float(x[2]), float(x[3]), int(x[4]))
            elif x[0] == 'steps':
                d['n'], d['wall'], d['mpi'] = int(x[1]), float(x[2]), float(x[3])
        ranks[r] = d

R = sorted(r for r in ranks if r not in IO)
n = ranks[R[0]]['n']
T_census = mean([ranks[r]['wall'] for r in R]) / n


def comp(d, ph):
    w, mpi, _ = d['phases'].get(ph, (0.0, 0.0, 0))
    return (w - mpi) / d['n']


# ---------------------------------------------------------------- byte model
def b_gs(s, f):          return 24 + s * (44 + 20 * f) + 8 * (s - 1)
def b_gs_stock(s, f):    return s * (24 + 44 + 20 * f)
def b_amul(f, r):        return 32 + r + 24 * f
def b_res(f, r):         return 48 + r + 24 * f
def b_fu(f):             return 72 + 15 * (108 + 20 * f) + 24 * 14
def b_res3(f, r):        return 128 + 3 * r + 24 * f


def per_rank(d, r):
    """Predicted change of compute time per step (s; negative = faster) for each idea, one rank."""
    L = d['levels']
    nC0, f0 = L[0]['nCells'], L[0]['nFaces'] / L[0]['nCells']
    ck, cw, cu = comp(d, 'solve:k'), comp(d, 'solve:omega'), comp(d, 'solve:U')
    kw = ck + cw
    share_gs8 = b_gs(8, f0) / (b_gs(8, f0) + 2 * b_res(f0, r))     # GS part of a k/omega solve
    share_fu = b_fu(f0) / (b_fu(f0) + 2 * b_res3(f0, r))          # sweep part of the U solve
    bw_kw = 2 * (b_gs(8, f0) + 2 * b_res(f0, r)) * nC0 / kw       # bytes/s of the scalar solves
    bw_fu = (b_fu(f0) + 2 * b_res3(f0, r)) * nC0 / cu             # bytes/s of the U solve
    am0 = L[0]['t_amul'] / n
    bw_am = b_amul(f0, r) * nC0 * L[0]['n_amul'] / L[0]['t_amul']  # bytes/s of the finest Amul
    lv = [(l['t_gs'] / n, l['t_amul'] / n, l['nFaces'] / l['nCells']) for l in L if l['nCells'] > 0]
    out = {}

    # compactAddr: Gauss-Seidel (GAMG levels, 2 sweeps; k/omega, 8 sweeps) and Amul (GAMG levels)
    out['compactAddr'] = -(
        sum(g * 2 * (3 + 2 * f) / b_gs(2, f) + a * (6 * f - 1) / b_amul(f, r) for g, a, f in lv)
        + kw * share_gs8 * 8 * (3 + 2 * f0) / b_gs(8, f0))

    # compact addressing inside the fused U solve (fusedVectorSolve v4)
    out['compactFusedU'] = -cu * share_fu * 15 * (3 + 2 * f0) / b_fu(f0)

    # folded per-sweep copy (stepReuse): GAMG levels (2 sweeps per call), k/omega (8 sweeps)
    out['gsFold'] = -(sum(g * 16 * (2 - 1) / b_gs(2, f) for g, a, f in lv)
                      + kw * share_gs8 * 16 * (8 - 1) / b_gs(8, f0))

    # fused Ux/Uy/Uz Gauss-Seidel v1 (shared coefficients/addressing, copy per sweep) vs 3 x 15
    # separate one-sweep calls
    out['fusedU_v1'] = -15 * (3 * (68 + 20 * f0) - (180 + 20 * f0)) * nC0 / bw_fu

    # fixed-iteration solvers (fvSolution renamed smoothSolver only for U, k, omega; GAMG stayed
    # stock, and the coarsest-level PCG keeps its default tolerance, so the census still shows 3
    # global sums per PCG iteration): per-iteration residual + norm removed (U 3 x 14, k 7,
    # omega 7 = 56 passes), with its global sum and its halo exchange
    out['fixedIter_bytes'] = -56 * b_res(f0, r) * nC0 / bw_am
    ar_lat = d['rest']['t_allreduce'] / d['rest']['n_allreduce']    # latency of a global sum
    halo = (L[0]['t_ifInit'] + L[0]['t_ifUpd']) / L[0]['n_ifInit']   # packing of a halo exchange
    out['fixedIter_sync'] = -56 * (ar_lat + halo)
    # cap: if every removed sum also removed the wait seen at a finest-level sum (imbalance)
    out['fixedIter_waitcap'] = -56 * L[0]['t_allreduce'] / L[0]['n_allreduce']

    # GaussSeidelNoCopy (fixed-iteration v1 era: every call one sweep, U unfused): no copy pass,
    # source read in the sweep, 1-byte flags per cell and per face: saves 15 - f bytes per sweep
    out['gsNoCopy'] = -(sum(g * 2 * (15 - f) / b_gs(2, f) for g, a, f in lv)
                        + (16 + 45) * (15 - f0) * nC0 / bw_kw)

    # cell-based Amul (lowerCSR) for the finest p products: one pass (r = 0), ownerStart and
    # losortStart, per face upper + upperAddr + lowerCSR + losort + lowerAddr, lowerCSR built
    # once per solve (25 products per step)
    b_csr = 40 + 28 * f0 + 28 * f0 / 25
    out['csrAmul'] = am0 * (b_csr - b_amul(f0, r)) / b_amul(f0, r)

    out['_bw'] = (bw_am, bw_kw, bw_fu,
                  b_gs(2, f0) * nC0 * L[0]['n_gs'] / L[0]['t_gs'] if L[0]['t_gs'] else 0)
    return out


lo = {r: per_rank(ranks[r], R_LO) for r in R}
hi = {r: per_rank(ranks[r], R_HI) for r in R}
compute = {r: (ranks[r]['wall'] - ranks[r]['mpi']) / n for r in R}
slow = sorted(R, key=lambda r: -compute[r])[:20]


def pred(key, sel=R):
    a = mean([lo[r][key] for r in sel]); b = mean([hi[r][key] for r in sel])
    return min(a, b), max(a, b)


# ---------------------------------------------------------------- sanity: achieved bandwidth
bw = [mean(lo[r]['_bw'][i] for r in R) / 1e9 for i in range(4)]
bwh = [mean(hi[r]['_bw'][i] for r in R) / 1e9 for i in range(4)]
print(f'census: {len(R)} ranks (I/O ranks excluded), {n} steps, {T_census:.4f} s/step with census on')
print('achieved bandwidth per rank, mean over ranks (GB/s; r = 0 .. 24 B/cell re-read):')
for i, name in enumerate(['finest Amul', 'k/omega solves', 'U solve', 'finest GAMG Gauss-Seidel']):
    print(f'  {name:26s} {bw[i]:5.2f} .. {bwh[i]:5.2f}   x104 ranks = '
          f'{104 * bw[i]:4.0f} .. {104 * bwh[i]:4.0f} GB/s per node (STREAM 492)')

# ---------------------------------------------------------------- predictions vs A/B
AB = [  # key, label, A/B, base s/step, measured change s/step (negative = faster), note
    ('fixedIter', 'Fixed-iteration solvers', 'fixediter/drivaer_4node_181329240', 1.9947, 1.6688 - 1.9947, ''),
    ('fusedU_v1', 'Fused Ux/Uy/Uz Gauss-Seidel', 'fusedu/drivaer_4node_181333833', 1.7739, 1.7070 - 1.7739, ''),
    ('gsFold', 'Folded Gauss-Seidel copy', 'stepreuse/bench_4node_181376795', 1.5694, 1.4993 - 1.5144,
     'difference of two bundle arms'),
    ('compactAddr', 'Compact addressing (GS + Amul)', 'compact/ab_4node_182180991', 1.26522, 1.24507 - 1.26522,
     'mean of 3 repeats'),
    ('compactFusedU', 'Compact addressing in fused U', 'compact/fusedu4_ab_4node_182181441', 1.24378,
     1.24188 - 1.24378, 'mean of 3 repeats'),
    ('gsNoCopy', 'GaussSeidelNoCopy (branchy)', 'gsnocopy/drivaer_1node_181331034', 6.9590, 7.1508 - 6.9590,
     '1-node A/B; prediction in % at 4n rates'),
    ('csrAmul', 'Cell-based Amul (gather)', 'csramul/bench_4node_181362805', 1.6121, 1.6513 - 1.6121, ''),
]
print()
print('| Idea | A/B | measured | predicted (all ranks) | predicted (20 slowest ranks) | verdict |')
print('|---|---|---:|---:|---:|---|')
rows = []
for key, label, ab, base, meas, note in AB:
    if key == 'fixedIter':
        a = pred('fixedIter_bytes'); s = pred('fixedIter_sync'); w = pred('fixedIter_waitcap')
        p = (a[0] + s[0] + w[0], a[1] + s[1]); ps = p
        extra = f' (bytes {100*a[1]/base:+.1f}..{100*a[0]/base:+.1f}%, syncs {100*s[0]/base:+.1f}%, waits at removed syncs up to {100*w[0]/base:+.1f}%)'
    else:
        p, ps = pred(key), pred(key, slow); extra = ''
    if key == 'gsNoCopy':     # prediction made at 4-node rates against the fixed-iteration-era step
        base_pred = 1.6688
    else:
        base_pred = base
    m = 100 * meas / base
    plo, phi = 100 * p[0] / base_pred, 100 * p[1] / base_pred
    slo, shi = 100 * ps[0] / base_pred, 100 * ps[1] / base_pred
    hit = (plo - 0.35 <= m <= phi + 0.35)
    sign_ok = (m < 0) == (phi < 0) or abs(m) < 0.35
    verdict = 'hit' if hit else ('direction right, size off' if sign_ok else 'MISS')
    rows.append((label, m, plo, phi, verdict))
    print(f'| {label} | `{ab}` | {m:+.2f}% | {plo:+.2f} .. {phi:+.2f}%{extra} | {slo:+.2f} .. {shi:+.2f}% | {verdict} |'
          + (f' {note}' if note else ''))


# ---------------------------------------------------------------- forward prediction
# Re-enabling the two function objects the Bonus arm defers (stock behaviour, Addendum F5).
# Stock per time step, full mesh (v2512 source):
#   wallShearStress::execute  devReff = (-nuEff) * dev(twoSymm(grad(U))), grad(U) from the cache:
#       nuEff = nut + nu (8 + 8 read, 8 + RFO 8 write) 32; unary minus 16 (reused tmp) .. 24;
#       twoSymm(gradU) 72 read + 48 write + RFO 48 = 168; dev 96 (reused) .. 144;
#       product 56 read + 48 write (reused) .. +48 RFO = 104 .. 152          -> 416 .. 520 B/cell
#   forceCoeffs::reset + calcForce/MomentCoeffs  forceCoeff == Zero, momentCoeff == Zero
#       (24 + RFO 24 each), forceCoeff = s*force(), momentCoeff = s*moment() (read 24, write
#       24 + RFO 24 into a tmp that is transferred)                              -> 240 B/cell
# Streaming field algebra run by all ranks at once: per-rank bandwidth between the node STREAM
# share (492/104 GB/s) and the finest Amul's measured rate.
STREAM_RANK = 492e9 / 104
fo = {}
for key, blo, bhi in (('wallShearStress', 416, 520), ('forceCoeffs fields', 240, 240)):
    t = []
    for r in R:
        nC0 = ranks[r]['levels'][0]['nCells']
        bw_fast = hi[r]['_bw'][0]
        t += [blo * nC0 / bw_fast, bhi * nC0 / STREAM_RANK]
    fo[key] = (mean(t[0::2]), mean(t[1::2]))
T_now = 1.2358    # current Bonus, 4 nodes x 150 steps (job 182182107)
print()
print('forward prediction (written before any measurement): cost of re-enabling the deferred')
print(f'function objects in the Bonus arm, 4 nodes, vs {T_now} s/step:')
tot = [0.0, 0.0]
for k, (a, b) in fo.items():
    tot[0] += a; tot[1] += b
    print(f'  {k:20s} +{1000 * a:5.1f} .. +{1000 * b:5.1f} ms/step   +{100 * a / T_now:4.2f} .. +{100 * b / T_now:4.2f}%')
print(f'  {"both":20s} +{1000 * tot[0]:5.1f} .. +{1000 * tot[1]:5.1f} ms/step   +{100 * tot[0] / T_now:4.2f} .. +{100 * tot[1] / T_now:4.2f}%')
