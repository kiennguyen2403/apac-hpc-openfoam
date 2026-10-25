# Performance models: simpleFoam occDrivAerStaticMesh 65M, 4 nodes x 104 ranks (team ya41)

Two models built from data we already had (no new runs), tested against the A/B results, plus one
forward prediction recorded before it is measured. Scripts: `tools/perfModel/bandwidth_model.py`,
`tools/perfModel/critical_path.py`; outputs: `bandwidth_model.txt`, `critical_path.txt`.

The scripts read the phase-census archive `diag/phase_census_4node_182140240.tar.gz` (per-rank timing from
`tools/gamgCensus`); job folders named below are on Gadi project storage and are not in this repository.
Usage: `python3 bandwidth_model.py [CENSUS_TAR]`, `python3 critical_path.py [CENSUS_TAR]`. A/B job script for the
forward prediction: `scripts/profile/compiler/fo_revert_ab_4node.sh`.

Input for both: the phase census `diag/phase_census_4node_182140240` (4 nodes x 104 ranks, 9 time
steps, per-rank self times of every kernel and phase with MPI and nested calls excluded, MPI time by
category). Build measured: Bonus tier just before compact addressing and the UCX setting
(1.3160 s/step with the census on, 1.2675 s/step without).

## 1. Memory-traffic model of the solver kernels

**Premise.** The time loop is memory-bandwidth bound: it streams 72% of the node's STREAM triad
bandwidth, and 72% of cycles have >= 4 reads outstanding (`vtune/membw_1node_182171696`). So a kernel's
time is proportional to the bytes it moves, and an optimisation that changes a kernel's bytes from B to
B' changes its time by t x (B' - B) / B, with t the kernel's measured time. **Nothing is fitted.**

**Bytes** are counted from the loops: 8-byte scalars, 4-byte labels; a write to a line not in cache
costs a read for ownership plus a write-back. f = faces per cell of the multigrid level (2.97 on the
finest level, up to ~5 on coarse levels). Per-face arrays are exact; whether per-cell arrays are
re-read inside one call is uncertain, so the model carries a range r = 0..24 B/cell.

| Kernel | Bytes per cell |
|---|---|
| Gauss-Seidel call of s sweeps, folded copy | 24 + s(44 + 20f) + 8(s - 1) |
| Gauss-Seidel, stock (copy every sweep) | s(24 + 44 + 20f) |
| Amul (matrix-vector product) | 32 + r + 24f |
| residual + its norm | 48 + r + 24f |
| fused Ux/Uy/Uz Gauss-Seidel, 15 sweeps | 72 + 15(108 + 20f) + 24 x 14 |
| compact addressing saves | s(3 + 2f) per Gauss-Seidel call, 6f - 1 per Amul |

**Achieved bandwidth per rank** (bytes / measured kernel time, mean over ranks): finest-level
Gauss-Seidel 6.6 GB/s, finest Amul 5.9-7.3, k/omega solves 4.7, U solve 3.6. x104 ranks this is
above the node's 492 GB/s STREAM triad for the multigrid kernels. Plausible reasons (not measured
separately): about a fifth of the rank-time is spent in MPI, so the ranks that are computing share the
bandwidth among fewer, and these read-dominated kernels beat the triad's 2:1 read/write mix.

**Predicted vs measured** (A/B results from same-job paired runs; noise about 0.35%):

| Idea | A/B | Measured | Predicted | Decision the model implies | Right? |
|---|---|---:|---:|---|---|
| Fixed-iteration solvers | `fixediter/drivaer_4node_181329240` | -16.3% | -9.4 .. -12.1% | do it | yes (size under-predicted) |
| Fused Ux/Uy/Uz Gauss-Seidel | `fusedu/drivaer_4node_181333833` | -3.8% | -5.1 .. -5.3% | do it | yes |
| Folded Gauss-Seidel copy | `stepreuse/bench_4node_181376795` (difference of two arms) | -1.0% | -1.0 .. -1.1% | do it | yes (size hit) |
| Compact addressing (GS + Amul) | `compact/ab_4node_182180991` (3 repeats) | -1.6% | -2.5 .. -2.7% | do it | yes |
| Compact addressing in fused U | `compact/fusedu4_ab_4node_182181441` (3 repeats) | -0.15% | -0.5% | marginal (at the noise floor) | yes (size within noise) |
| Cell-based Amul (gather loop) | `csramul/bench_4node_181362805` | +2.4% | 0 .. +1.0% | skip | yes |
| GaussSeidelNoCopy (per-face branches) | `gsnocopy/drivaer_1node_181331034` (1 node) | +2.8% | -2.2% | do it | **no** |

- **6 of 7 decisions right.** The miss is instructive: NoCopy saves the same copy pass the later
  fold saves, but with a branch per face; the branch-free fold won (-1.0%), the branchy version lost.
  Bytes are necessary, not sufficient: a kernel that gains branches stops being bandwidth bound.
- **Size.** For pure layout changes the model over-predicts by 1.1-1.7x (realisation 60-90% of the
  byte savings): the sequential streams shrink, the irregular neighbour accesses do not.
- **Fixed-iteration solvers** are under-predicted: bytes explain 8.7-9.0%, the 56 removed global sums
  and halo exchanges 0.7%, waits at those sums up to 2.5%; the remaining ~4% is not explained (most
  likely more waiting at the removed global syncs than the finest-level sums suggest).
- **Caveat.** These A/B results were known when the byte counts were written, so this is a
  retrospective test, even with nothing fitted. Section 3 is the real test.

## 2. Critical-path decomposition of a time step

Each step is a chain of segments (each census phase; the pressure solve split into its multigrid
levels), each ending in a synchronisation. To first order step = sum of segment maxima + communication:

| Part | ms/step | share |
|---|---:|---:|
| Average compute (mean over ranks) | 1040.6 | **79.1%** |
| Whole-step imbalance (slowest rank - mean) | 66.2 | 5.0% |
| Per-segment imbalance (different ranks slowest in different segments; lower bound) | 61.8 | 4.7% |
| Communication + imbalance finer than a segment (upper bound) | 147.4 | 11.2% |

Largest per-segment imbalances: finest multigrid level 14.6 ms, grad(U) 11.3, second level 9.4, U solve
8.6, laplacian(U) 7.4, omega wall function 7.0 ms/step. MPI time per step: 277 ms on average, 225 ms
on the 20 slowest ranks (Allreduce 103, Wait 74, posting 31, other collectives 16).

**Ceilings and what we found:**
- **Perfect balance in every segment: at most 9.7%.** We searched 82 scotch/hierarchical partitions,
  KaHIP configurations and weighted partitions: nothing beat scotch (4 104)/(1 0.5). Half of the
  ceiling is per-segment imbalance, which a single-constraint partitioner (scotch, KaHIP) cannot target:
  it balances total cells, not each multigrid level or phase separately.
- **Zero-cost communication: at most 11.2%** (includes imbalance finer than a segment). Transport,
  protocol and collective tuning measured as noise; the one communication-side gain was the UCX
  setting (-2.7%), which removed a compute slowdown on ranks with off-node neighbours.
- **79% is bandwidth-bound compute**, so the remaining lever was bytes (Section 1), which is where the
  later gains came from.

## 3. Forward prediction (recorded 2026-10-25T12:00:22Z, before any measurement)

The compliant Bonus arm re-enables two function objects that the current Bonus arm defers
(Addendum F5: function objects are part of the work). Stock per-step work, full mesh (v2512 source):
`wallShearStress` builds devReff = (-nuEff) dev(twoSymm(grad(U))) (416-520 B/cell); `forceCoeffs`
zeroes and refills the forceCoeff and momentCoeff fields (240 B/cell). Streaming field algebra, all
ranks at once: per-rank bandwidth between the STREAM share (4.7 GB/s) and the finest Amul's rate.

| Re-enabled | Predicted cost, 4 nodes (vs 1.2358 s/step) |
|---|---:|
| wallShearStress | +9.0 .. +17.3 ms/step (+0.73 .. +1.40%) |
| forceCoeffs fields | +5.2 .. +8.0 ms/step (+0.42 .. +0.64%) |
| both | **+14.2 .. +25.2 ms/step (+1.15 .. +2.04%)** |

Test: same-job A/B at 4 nodes, current Bonus vs the same with `STEP_REUSE_OFF=wss,fco` (and a middle
arm with `wss` only), fixed arm order, 3 repeats. Checksums of the model and its output at the time of
writing: `prediction_written_utc.txt`.

**Result** (job 182232191, measured 2026-10-25 12:15-12:18 UTC, after the prediction was recorded; one
repeat per arm at the user's request, same nodes, fixed arm order, p residuals bit-identical in all arms;
folder `fo_revert_ab_4node_182232191`):

| Arm | s/step (20 steps) | Measured vs cur | Predicted | Consistent (+-0.35% noise)? |
|---|---:|---:|---:|---|
| cur (both deferred) | 1.2305 | - | - | |
| wss (stock wallShearStress) | 1.2480 | **+1.43%** | +0.73 .. +1.40% | yes (top of the range) |
| wssfco (+ stock forceCoeffs fields) | 1.2579 | **+2.23%** | +1.15 .. +2.04% | yes (0.19 above the range) |
| forceCoeffs fields alone (wssfco - wss) | | +0.80% | +0.42 .. +0.64% | yes (0.16 above the range) |

The forward prediction holds within the noise of a single repeat, with all three results at the upper end.
The upper end is the bound that assumes every rank streams at its STREAM share at the same time, which
is how function objects run: all ranks together after the turbulence update, with no MPI waits to free
bandwidth. Unlike the layout changes in Section 1 (60-90% realised), plain streaming field algebra
converts its bytes to time fully. Cost of the compliant Bonus arm from this item: about +2.2%.
