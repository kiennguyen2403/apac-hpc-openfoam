# OpenFOAM v2512 performance work — APAC HPC-AI 2026 (occDrivAerStaticMesh 65M)

Tools, solver code and job scripts from tuning `simpleFoam` on the 65M-cell
occDrivAerStaticMesh case on NCI Gadi (Sapphire Rapids, `normalsr`), OpenFOAM v2512
built with Intel oneAPI ICX 2025.3 and Intel MPI 2021.17.

## Results (4 nodes × 104 ranks, 150 steps, best of 3, arms interleaved on the same nodes)

| Setup | s/step | vs Basic | Cd at step 150 |
|---|---:|---:|---|
| Basic: ICX-SPR + Intel MPI (tuned Allreduce, scatter pinning), scotch `(4 104)/(1 0.5)`, glibc malloc tuning | 1.9765 | – | 0.28738062413153675 |
| Bonus: Basic + `fixedIterSolvers` | 1.6684 | −15.6% | identical |
| Bonus: + `fusedVectorSolve` v1 | 1.6066 | −18.7% | identical |
| Bonus: + `fusedVectorSolve` v2 | 1.5919 | −19.5% | identical |
| Bonus: + `stepReuse`, `fixedIterSolvers` v2, `fusedVectorSolve` v3, cache `grad(p)`, `grad(k)`, `grad(omega)` | 1.4815 | −25.0% | identical |
| Bonus: + `stepReuse` divDevReff and kOmegaSST per-cell fusion | 1.4242 | −27.9% | identical |
| Bonus: + `stepReuse` laplacian, convection, gradient fusion, lazy forceCoeffs fields | 1.3060 | −33.9% | identical |
| Bonus: + `simpleFoamFused` | 1.2675 | −35.9% | identical |
| Bonus: + fused cellLimited gradient, UCX zero-copy/rendezvous off | 1.2580 | −36.4% | identical |
| **Bonus: + `compactAddr`** | **1.2358** | **−37.5%** | identical |

Basic and the first Bonus row ran in the same job (181332282); the v1 and v2 rows are separate jobs
(181334343, 181363792, 181376794, 181390148, 181413566, 182135721, 182166454, 182182107), each checked against the
same Cd. The Basic row predates the UCX setting (both arms of `scripts/final/final_runs.sh` now use it). Both solver changes are bit-identical to stock OpenFOAM: every solver log line, field and force
coefficient matches.

2 and 1 nodes (current setup, Basic and Bonus in the same job, Cd bit-identical):

| Configuration | Basic s/step | Bonus s/step | vs Basic | Cd |
|---|---:|---:|---:|---|
| 2 nodes × 104, 80 steps (job 182181634) | 4.0081 | 2.5402 | −36.6% | 0.33772590153571347 |
| 1 node × 104, 40 steps (job 182181635) | 8.1711 | 5.2129 | −36.2% | 0.43805628753637488 |

**Contest addendum (2026-09-25) compliance.** The addendum keeps `fvSolution.fixedIter`, `fvSchemes` and
`controlDict.noWrite` (including its function objects) as shipped in every tier, and limits the Basic case edits
to `nCores`, `decompositionMethod`, `nHierarchical`, the method's coefficients and `endTime`. The numbers above
used the team's `controlDict` (`runTimeModifiable no`), and the Bonus arm also renamed solvers and extended the
cache in `fvSolution` and deferred two function objects. Status:

- `scripts/final/final_runs.sh` now builds `controlDict` from the shipped `controlDict.noWrite` (only `endTime`
  changed) and stops if `fvSolution.fixedIter`, `fvSchemes` or `controlDict.noWrite` differ from the shipped case,
  or if `caseDefinition` differs outside `nCores` / `decompositionMethod` / `nHierarchical` (the scotch
  `domains` / `domainWeights` now sit in `decomposeParDict`'s `scotchCoeffs`).
- **Compliant Basic, 4 nodes × 150 steps: 1.9341 s/step** (job 182187466, one run, with the UCX setting),
  Cd 0.28738062413153675, identical to the earlier Basic.
- Re-enabling the deferred `wallShearStress` and `forceCoeffs` fields costs +2.23% in the Bonus arm (job
  182232191). Moving the fixed-iteration and cache switches from `fvSolution` into the preloaded libraries, and
  the compliant Bonus finals, are still to do.

## `tools/fixedIterSolvers` — solver plugin (main result)

With fixed-iteration controls (`tolerance 0; relTol 0;` plus `maxIter`), convergence can
never be reached, yet the stock solvers still test for it every iteration:

- `smoothSolver` computes a full residual after every sweep (a new field, one extra matrix
  pass and a global reduction);
- `PCG` computes a residual norm (a global reduction) every iteration.

`smoothSolverFixed` and `PCGFixed` do the same sweeps and iterations and compute these once,
after the last one. When convergence is possible they behave exactly like the originals.

Results are **bit-identical** to the stock solvers (every solver log line, every field and
force coefficient), measured paired on the same nodes, 10 steps:

| Run | Stock | Plugin | Change |
|---|---:|---:|---:|
| 1 node × 104 ranks | 8.165 s/step | 6.988 s/step | −14.4% |
| 4 nodes × 104 ranks | 1.995 s/step | 1.669 s/step | −16.3% |

Build with the same flags as the OpenFOAM build (required for bit-identity, since `PCGFixed`
carries its own copy of PCG's arithmetic loops):

```bash
FOAM_EXTRA_CXXFLAGS="-xSAPPHIRERAPIDS -fp-model=precise" wmake libso
```

v2: with a GaussSeidel smoother, `smoothSolverFixed` does all the sweeps of its iteration loop in
one smoother call (n sweeps in one call are exactly n calls of one sweep), so the smoother can fold
its per-sweep copy (see `tools/stepReuse`).

`GaussSeidelNoCopy` (same library) is a GaussSeidel smoother that avoids the per-sweep
`bPrime = source` copy. It is bit-identical but measured **2.8% slower** (the per-face branch
costs more than the streaming copy), so it is not used.

Enable in a case:

```
// system/controlDict
libs ("<path>/libfixedIterSolvers.so");

// system/fvSolution: smoothSolver -> smoothSolverFixed, coarsest-level PCG -> PCGFixed
```

## `tools/fusedVectorSolve` — one GaussSeidel sweep for Ux, Uy, Uz

OpenFOAM solves the three velocity components one after another. They are independent and share
the matrix off-diagonals and addressing; only the boundary-corrected diagonal (different per
component next to `slip`/`symmetry` patches), the source and the interface coefficients differ.
This library replaces `fvMatrix<vector>::solveSegregated`: each component is set up exactly as in
OpenFOAM (own diagonal, source, interface correction, initial residual, norm factor), then one
GaussSeidel sweep updates all three components, reading the shared matrix once instead of three
times. It only engages for `smoothSolver`/`smoothSolverFixed` + `GaussSeidel` with
`tolerance 0; relTol 0;` on a 3D mesh; anything else runs a verbatim copy of the original code.

`libfiniteVolume` exports the generic instantiation as a weak symbol called through the PLT, so the
library is activated with `LD_PRELOAD` (no OpenFOAM rebuild):

```bash
mpirun ... -genv LD_PRELOAD <path>/libfusedVectorSolve.so simpleFoam -parallel
# FUSED_VECTOR_SOLVE_REPORT=1 prints (stderr, once) whether the fused or generic path ran
```

Paired, same nodes, 10 steps, on top of `fixedIterSolvers`: −2.77% (1 node), −3.77% (4 nodes),
bit-identical. Build with the same `FOAM_EXTRA_CXXFLAGS` as above.

v2 also fuses the passes around the sweeps: the three components' A·ψ and row sums (`Amul`, `sumA`,
for the initial residual and norm factor) are computed in one pass over the faces, and so are the
three final residuals, instead of nine separate matrix passes. Global reductions stay per component,
so the MPI summation order is unchanged. Paired against v1 (4 nodes, 20 steps): **−1.22%**,
bit-identical. Measured in a 150-step final run: 1.5919 s/step (v1: 1.6066).

v3 runs all sweeps as one sequence and folds the per-sweep `bPrime = source` copy into the previous
sweep, as `tools/stepReuse` does for the scalar GaussSeidel smoother.

`test/` holds the bit-identity check used for every change: `make_box_case.sh` builds a 3D box case
with the same boundary types as occDrivAerStaticMesh (slip, symmetry, wall functions, kOmegaSST), and
`identity_test.sh` compares stock OpenFOAM with the preloaded library, serial and on 4 ranks
(`qsub -v WORK=<dir> identity_test.sh`). Options: `SOLVER_SED`/`LIBS` (solver plugins in the test
arm), `CACHE_KOMEGA=1`, `WSS=1` (wallShearStress function object, output compared), `DEBUG_CACHE=1`
(log solution cache use) and `BOX_SKEW=1` (sheared, non-orthogonal box, so the corrected laplacian
and snGrad terms compute their gradient correction).

## `tools/stepReuse` — work done twice per step, done once

One `LD_PRELOAD` library with three replacements, each bit-identical to stock OpenFOAM:

- **GaussSeidel**, `GaussSeidelSmoother::smooth`: every sweep starts with `bPrime = source`, a
  separate streaming pass. The copy for the next sweep now happens inside the current one, right after
  `bPrime[celli]` is read for the last time (later cells only write `bPrime[j]` for `j > celli`). Same
  values, same operations, one memory pass less per sweep. Applies to p (GAMG smoothing), and with
  `smoothSolverFixed` v2 to k and omega; `fusedVectorSolve` v3 does the same for U.
- **kOmegaSST**, `kOmegaSSTBase<...>::correct` (explicit specialisation): `fvc::surfaceIntegrate(phi)`
  was computed three times (for `divU` and in the `bounded` correction of both convection terms) and
  the upwind convection coefficients of omega and k, identical because the weights are `pos0(phi)`,
  twice. Each is now built once; only the boundary coefficients are built per field.
- **kOmegaSST per-cell maths**: S2, GbyNu0, CDkOmega, F1, F23 (twice per step), GbyNu, gamma, beta,
  DkEff, DomegaEff and nut were each built from one field operation (and one temporary field) per
  arithmetic operation, about 60 cell-sized temporaries per step. Each is now one loop over the cells
  and one per patch. Every value goes through the same scalar operations in the same order, with FMA
  contraction off in these loops (`#pragma clang fp contract(off)`, checked to be honoured by ICX; the
  stock field operations are one operation per loop, so nothing is contracted there), and the tensor
  functions are the same inline functions. Plain `kOmegaSST` with `F3 off` only.
- **divDevReff**, `linearViscousStress<...>::divDevRhoReff`: the explicit term
  `fvc::div(nuEff*dev2(T(grad(U))))` made three tensor temporaries, a face field and two passes, and
  `nuEff()` was evaluated twice. Now `X = nuEff*dev2(T(gradU))` is one cell loop and the interpolation,
  `Sf &` and surface integral one face loop, adding to owner and neighbour in face order and then the
  patches, as `surfaceIntegrate` does, so every cell sum is unchanged.
- **Laplacian** (`laplacianScheme<Type, scalar>::fvmLaplacian` with a volume gamma: U, k, omega and p):
  Gauss linear corrected made the interpolated gamma, gamma*|Sf|, the correction and gamma*|Sf|*correction
  as face fields, then a surface integral. Now one face loop does the interpolation, coefficients, diagonal
  (in negSumDiag order) and correction flux (in surfaceIntegrate order); the interpolation expression is the
  stock one, left to the same contraction. The flux correction p needs is stored as before.
- **U convection** (`gaussConvectionScheme<vector>::fvmDiv` with `linearUpwindV`): weights, coefficients and
  the limited correction (the linearUpwindV code) in one face loop, added to the cells directly instead of a
  correction face field, a flux face field and a surface integral.
- **Gauss gradient** (`gaussGrad<Type>::calcGrad`, both the returning and the in-place version used by the
  solution cache): the interpolation inside the gradient's face loop.
- **forceCoeffs**: the two full-mesh coefficient fields are only computed with `writeFields` on; the written
  coefficients come from the summed patch forces (`execute()` is replaced too, since the library inlines
  `reset()` into it).
- **cellLimited gradient of U** (`cellLimitedGrad<vector, minmod>::calcGrad`): the neighbour max/min is taken
  in the Gauss gradient's face loop, max/min/limiter are initialised in one loop and the deltas max - U,
  min - U are formed where the limiter uses them (stock: two copies, a separate max/min face loop, two
  subtraction passes). Max, min and the limiter are min/max selections, exact in any order. −0.37% at 4
  nodes (job 182137125), within noise; bit-identical.
- **wallShearStress** function object: evaluated `devReff()` over the whole mesh every step but only
  writes at write time (never, in a benchmark run). `execute()` now defers to `write()`, which runs the
  stock `execute()` and `write()` (through `dlsym(RTLD_NEXT)`).

Switch parts off with `STEP_REUSE_OFF=gs,komega,kfuse,divdev,lap,conv,grad,climit,fco,wss`;
`STEP_REUSE_REPORT=1` prints which path ran.

Also: `grad(k)` and `grad(omega)` in the fvSolution `cache` list. They are computed for `CDkOmega` and
again by the corrected laplacian terms; the stock solution cache reuses them (no code needed), like
`grad(p)`, which uncached is computed five times per step (box test with `DEBUG_CACHE=1`).

Paired, 4 nodes, 20 steps, same nodes (job 181376137), bit-identical:

| Arm | s/step | vs base |
|---|---:|---:|
| base: Bonus with `fusedVectorSolve` v2 and `grad(p)` cache | 1.6202 | – |
| + `cache { grad(k); grad(omega); }` | 1.6098 | −0.64% |
| + `stepReuse`, `fixedIterSolvers` v2, `fusedVectorSolve` v3 | 1.5503 | −4.31% |

`grad(p)` in the cache list alone: −2.30% (job 181366105, bit-identical).

Attribution (job 181376795, 4 nodes, 20 steps, same nodes): all −4.47%; without the kOmegaSST
part −2.97% (so ~1.5% from it); without the GaussSeidel fold −3.51% (~1.0%); the cache entries
~0.6% (above). Measured over 150 steps on top of `grad(p)` caching: 1.4815 s/step (job 181376794).

The divDevReff and kOmegaSST per-cell fusions, on top of that (job 181379590, 4 nodes, 20 steps, same
nodes, bit-identical):

| Arm | s/step | vs base |
|---|---:|---:|
| base: the setup measured at 1.4815 | 1.4959 | – |
| + divDevReff fused | 1.4744 | −1.43% |
| + kOmegaSST per-cell maths fused | 1.4431 | −3.53% |

The laplacian, convection, gradient and forceCoeffs changes, on top of the setup measured at 1.4242 (job
181412786, 4 nodes, 20 steps, same nodes, bit-identical):

| Arm | s/step | vs base |
|---|---:|---:|
| base: the setup measured at 1.4242 | 1.4382 | – |
| + convection, gradient, forceCoeffs | 1.3723 | −4.58% |
| + laplacian | 1.3119 | −8.78% |

## `tools/simpleFoamFused` — simpleFoam with fused momentum/pressure steps

A modified `simpleFoam` (its momentum and pressure steps live in the solver executable, so they cannot be
replaced by preloading). Bit-identical to `simpleFoam` (box identity test, serial and 4 ranks: every field,
solver line, wall shear stress and force coefficient; job 182135541). Paired, 4 nodes, 20 steps, same nodes,
with the production libraries in both arms (job 182135600): `simpleFoam` 1.3453 s/step, `simpleFoamFused`
1.3087 s/step (**−2.72%**), solver lines and force coefficients identical.

- `solve(UEqn == -fvc::grad(p))` deep-copies the momentum matrix to add `V*(-grad(p))` to the source. Here
  that source is built separately and swapped into `UEqn` for the solve, then swapped back (the solve leaves
  the matrix unchanged).
- `phiHbyA = fvc::flux(HbyA)` plus the SIMPLEC term `interpolate(rAtU - rAU)*snGrad(p)*magSf`: one face loop
  instead of the flux, interpolation, snGrad, correction and product face fields.
- `phi = phiHbyA - pEqn.flux()`: one face loop instead of the flux and difference face fields.
- `HbyA -= (rAU - rAtU)*grad(p)` and `U = HbyA - rAtU*grad(p)`: one cell loop each.

The two expressions the stock code leaves to FMA contraction (the linear interpolation and the matrix face
flux) are written the same way; every other operation is a separate statement, as it is a separate field
operation in the stock code. Falls back to the stock code (`stock/` holds the v2512 originals) unless the case
is SIMPLEC, has no MRF zones, has a fixed-value pressure patch, and uses linear interpolation and corrected
snGrad; `SIMPLEFOAM_FUSED_OFF=1` forces the stock code. Prepared: `identity_test.sh APP=<executable>` and
`scripts/profile/compiler/sff_4node.sh` (4 nodes, `simpleFoam` vs `simpleFoamFused`). The bonus arm of
`scripts/final/final_runs.sh` runs `simpleFoamFused`.

## `tools/compactAddr` — compact solver addressing

The time loop is memory-bandwidth bound (see *Memory roofline*), so the Gauss-Seidel smoother and `Amul`
read less addressing: `compactAddr` (LD_PRELOAD, **first** in the list) replaces
`GaussSeidelSmoother::smooth` (keeping `stepReuse`'s folded source copy) and `lduMatrix::Amul`. Per ldu level
(cached by `lduAddressing` pointer, `compactAddressing.H`) it stores the neighbour as a `uint16` offset from
the owner and the faces per cell as `uint8`, instead of `int32` `upperAddr`/`ownerStartAddr`; `Amul` no longer
reads `lowerAddr` (faces are owner-ordered, so the owner is the cell whose face range is walked). Faces are
visited in increasing index with the same statements, so results are bit-identical (box identity test, serial and
4 ranks, full production setup: job 182180957). Levels that do not fit fall back to the next definition;
`STEP_REUSE_OFF=compact` disables it.

Paired, 4 nodes, 20 steps, 3 repeats in fixed order (job 182180991): 1.2652 → 1.2451 s/step (**−1.6%**, every
compact run faster than every baseline run, pressure residuals identical). The same compaction inside the fused U
solve (`tools/fusedVectorSolve/v4build`, job 182181441) gave −0.15% (noise): the fused solve already reads the
addressing once for three components. Production keeps `fusedVectorSolve` v3.

## Memory roofline (`tools/memBandwidth`)

`scripts/profile/compiler/membw_1node.sh` (1 node, steps 2–12, `perf stat -I 1000` user-mode counters on every
rank; VTune's microarchitecture analysis needs a sampling driver that Gadi does not load):

- STREAM triad on all 104 cores (`stream_triad.c`): **492 GB/s** per node.
- The solver streams **353 GB/s** into L2 (72% of STREAM, MPI time included; 86% of the time at 300–425 GB/s),
  1.9 TB per step = 29 KB per cell per step.
- 71.8% of cycles have ≥4 data reads outstanding (bandwidth-limited), 5.2% have 1–3 (latency-limited).
- DTLB page walks: 1.2% of cycles, so huge pages cannot help.
- Store traffic (`rfo_1node.sh`, `analyse_rfo.py`): demand read-for-ownership misses are 9.4% of the L2 fills and
  demand RFOs sent beyond L2 12.1% (lower bounds of what streaming stores could save); dirty write-backs equal
  99% of the fills (upper bound).

## Performance models (`tools/perfModel`)

Built from the per-rank census (no new runs); details and outputs in `tools/perfModel/README.md`.

- `bandwidth_model.py`: kernel time proportional to the bytes it moves, bytes counted from the loops, nothing
  fitted. Against seven A/B results it gets 6 of 7 decisions right; the miss is `GaussSeidelNoCopy`, which saves
  the same bytes as the later folded copy but adds a branch per face. It over-predicts pure layout changes by
  1.1–1.7× and explains 60–75% of the fixed-iteration gain.
- `critical_path.py`: a step is 79% average (bandwidth-bound) compute, 5.0% whole-step imbalance, ≥4.7%
  per-segment imbalance and ≤11.2% communication, so perfect partitioning could save at most 9.7% and
  free communication at most 11.2%.
- Forward prediction, recorded before the run: re-enabling the deferred function objects costs +1.15 .. +2.04%;
  measured +2.23% (job 182232191, one repeat), consistent within noise.

## `tools/decompScore` — decomposition screening

- `decompScore` (OpenFOAM utility): reads the mesh once, then forks one child per candidate
  `decomposeParDict` and writes per-rank statistics (cells, processor faces, boundary faces,
  neighbours, off-node faces) without writing processor directories. Verified identical to
  `decomposePar`. Forking matters: running several scotch decompositions in one process
  does not reproduce `decomposePar`.
- `scotchStrat`: builds valid scotch strategy strings from scotch's own flags.
- `score_decomp.py`, `refit_validate.py`: per-rank cost model fitted on IPM data; ranks
  candidates by their predicted slowest rank.

Finding: per-step cost per rank tracks neighbour count; no scotch setting beat
`domains (4 104); domainWeights (1 0.5);` for this case. Note OpenFOAM truncates
`domainWeights` to integer ratios (`int(w/min)`), so `(1 0.4)` equals `(1 0.5)`.

## Load-imbalance diagnosis

About 19% of each step is ranks waiting for the slowest one. `scripts/profile/compiler/diag_imbalance_4node.sh`
measures per-rank compute (APS elapsed − MPI, 25 minus 5 steps) twice: with the normal rank placement and
with node order and core order reversed. `tools/decompScore/analyse_imbalance.py` compares the two.

- Per-rank time follows the **subdomain** (correlation +0.92 across placements), not the hardware (no
  node, socket or NUMA effect). The slowest rank is ~7% above the median.
- One scotch node group (ranks 0–103) is ~11% cheaper per rank than the others, with the same cell count.
- `decompScore` now also reports internal faces, wall faces and wall cells per rank;
  `fit_cell_cost.py` fits per-rank cost on them: cost follows **faces** (~7.5 µs per face per step) plus
  ~6.5 µs extra per wall face; cells add nothing once faces are counted (R² 0.39).
- `tools/cellWeights` writes integer per-cell weights for `decomposePar`'s `weightField` (faces-based, or
  scaled per rank by measured cost). OpenFOAM converts weights with `int(w/min)`, so the tool sets one
  cell to 1 and dithers the rest; otherwise most weights would be truncated away.
- `tools/csrAmul` (LD_PRELOAD) builds `lowerCSR` for the pressure matrix so `lduMatrix::Amul` takes its
  cell-based (gather) loop. `losort` is stable, so this is bit-identical.
- `tools/smtprobe` + `smt_probe_1node.sh`: checks whether a second rank per core (hyperthread) can be used.

Root cause, found with `tools/gamgCensus` (LD_PRELOAD timing wrappers, per rank: GAMG levels, step phases,
MPI time, neighbours, rusage; `analyse_*.py`):

- GAMG coarse levels are balanced in size; their face counts and the processor interfaces vary, and the pressure
  solve explains only 28–39% of the per-rank variance. A common per-rank efficiency factor slows every loop on
  some ranks. Not locality (CuthillMcKee band ~10k cells), not denormals (FTZ/DAZ bit-identical, zero FP assists,
  job 182140425).
- Extra L3 misses track off-node processor faces. Shifting the rank-to-node mapping by 52 (job 182165344)
  moved the slowdown with the off-node faces: ~5 µs per off-node face per step, inside loops without MPI calls.
- UCX registered the solver's buffers for zero-copy/rendezvous transfers to other nodes.
  `UCX_ZCOPY_THRESH=inf UCX_RNDV_THRESH=inf` (eager copies through UCX's pre-registered buffers) removed the
  penalty (job 182165954: MPI 257 → 213 ms/step, voluntary context switches 11 → 0.3 per step) and gave
  **−2.7%** within one job (3 repeats, job 182166167), bit-identical. Both arms of `final_runs.sh` use it.

## Tried and rejected

| Change | Result |
|---|---|
| `GaussSeidelNoCopy` smoother (no per-sweep source copy) | bit-identical, 2.8% slower |
| GAMG `scaleCorrection false` for p | 6.7% faster, but Cd −0.608% at step 150 (Bonus limit ±0.5%) |
| Single-reduction coarsest-level PCG | not built: APS showed `Allreduce` is ~89% load-imbalance wait, ~0.4% to gain |
| KaHIP (team partition, max 37 neighbours, half the inter-node faces) | −0.36% vs scotch `(1 0.5)`, within noise |
| scotch domainWeights / strategies / hierarchical (82 candidates, per-rank cost model) | nothing better than `(4 104)/(1 0.5)` |
| Weighted decomposition (`weightField`, faces-based / measured) | +5.6% / +8.3%: weights worsen scotch's cuts (max processor faces 24k → 40k); the cheap group's advantage is not a per-cell property |
| Cell-based `Amul` for p (`csrAmul`, cell loops suggested by Culpo, PRACE 2012) | bit-identical, 2.4% slower: the extra index lookup and per-step CSR rebuild cost more than the scattered writes it avoids (CuthillMcKee ordering keeps those writes in cache); a cell-based GaussSeidel was therefore not built |
| `GAMGFused` (`tools/gamgFused`): copy of GAMG with the residual, restriction and finest-level scaling fused | bit-identical, +0.2% (noise): the fused cell loops are cheap next to the face loops |
| Hybrid MPI + OpenMP | not built: all 104 cores per node already run a rank, v2512 has no OpenMP in the simpleFoam path, the GaussSeidel smoother is sequential; ceiling about the MPI share (~10%) |
| Direct assembly of the omega and k equations (one face loop for convection, laplacian and correction, one cell loop for diagonal and source, instead of ~9 matrices combined pairwise) | bit-identical, +1.16% at 4 nodes (job 182136701): the face loop's scatter into three extra cell accumulators costs more than the streaming pairwise matrix operations it replaces |
| Cp (`pressure` function object) evaluated at write only | not done: the patch probes sample `cp`, so deferring it is only correct because of the function-object order; ~0.15% |
| Two GaussSeidel sweeps per pass over memory (second sweep trailing the first through a cache window) | not built: keeping the results exact needs processor-boundary cells last and a window of about two boundary layers; from the per-rank processor faces that is 2.2–4.2 MB for p and 4.5–8.8 MB for the fused U sweep per rank, larger than the 2 MB L2, and largest on the slowest ranks. A version that fits (current ordering, stale halo values in every second sweep) changes the smoother |
| Two ranks per core (hyperthreading) | not available: a whole-node `normalsr` job's cpuset is `0-103`, one hyperthread per core |
| GAMG processor agglomeration (masterCoarsest 104/26, procFaces 50/400), current setup | +0.80% / −0.24% / −0.49% / +0.47%: noise, and not bit-identical |
| Weighted scotch v2 (`decompScore` candidates `name dict A B C`, w = A·internal + B·wall + C) | every candidate predicted 2.6–8.1% slower: weights inflate max processor faces 24k → 33–48k |
| KaHIP tuning (strong/eco imbalance 0.3% / 0.1%, cell weights, flat 416-way, fast) scored by `decompScore` + census cost model | best variants equal scotch on the slowest rank; the team partition measured −0.4..−0.9% but moves Cd −0.45% at step 20 |
| UCX transports (RC, DC, legacy protocols) and HCOLL/SHARP | within ±0.6% (noise); HCOLL/SHARP did not engage |
| Compact pinning / rank reordering to keep halos within a socket | with the UCX setting, cross-socket faces cost the same as same-socket faces (census fit) |
| Huge pages for the heap | page walks are 1.2% of cycles; jemalloc with THP was −2.5% on one node, +92% on another |
| Compact addressing in the fused U solve (`fusedVectorSolve` v4) | bit-identical, −0.15% (noise) |
| `stepReuse` deferral of `wallShearStress` to write time and of the `forceCoeffs` coefficient fields | output identical, but the work is skipped; dropped for addendum compliance (function objects are part of the defined work), costs +2.23% |

## `scripts/` and `wmake-rules/`

PBS job scripts (builds, profiling with perf / VTune / APS / IPM, benchmarks) and the extra
wmake rule files (`OptSPR` and variants) used with the OpenFOAM tree. Paths are specific to
Gadi project storage.

## License

GPL-3.0, as OpenFOAM: `PCGFixed.C` is a modified copy of OpenFOAM's `PCG.C`. See `LICENSE`.
