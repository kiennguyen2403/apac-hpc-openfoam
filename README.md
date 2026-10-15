# OpenFOAM v2512 performance work — APAC HPC-AI 2026 (occDrivAerStaticMesh 65M)

Tools, solver code and job scripts from tuning `simpleFoam` on the 65M-cell
occDrivAerStaticMesh case on NCI Gadi (Sapphire Rapids, `normalsr`), OpenFOAM v2512
built with Intel oneAPI ICX 2025.3 and Intel MPI 2021.17.

## Results (4 nodes × 104 ranks, 150 steps, best of 3, arms interleaved on the same nodes)

| Setup | s/step | vs Basic | Cd at step 150 |
|---|---:|---:|---|
| Basic: ICX-SPR + Intel MPI (tuned Allreduce, scatter pinning), scotch `(4 104)/(1 0.5)`, glibc malloc tuning | 1.9765 | – | 0.28738062413153675 |
| Bonus: Basic + `fixedIterSolvers` | 1.6684 | −15.6% | identical |
| Bonus: + `fusedVectorSolve` | 1.6066 | −18.7% | identical |

Basic and the first Bonus row ran in the same job (181332282); the last row is a separate job (181334343),
checked against the same Cd. Both solver changes are bit-identical to stock OpenFOAM: every solver log line, field and force
coefficient matches.

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
bit-identical. The 1.6066 s/step final result above was measured with v1.

`test/` holds the bit-identity check used for every change: `make_box_case.sh` builds a 3D box case
with the same boundary types as occDrivAerStaticMesh (slip, symmetry, wall functions, kOmegaSST), and
`identity_test.sh` compares stock OpenFOAM with the preloaded library, serial and on 4 ranks
(`qsub -v WORK=<dir> identity_test.sh`).

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

## Tried and rejected

| Change | Result |
|---|---|
| `GaussSeidelNoCopy` smoother (no per-sweep source copy) | bit-identical, 2.8% slower |
| GAMG `scaleCorrection false` for p | 6.7% faster, but Cd −0.608% at step 150 (Bonus limit ±0.5%) |
| Single-reduction coarsest-level PCG | not built: APS showed `Allreduce` is ~89% load-imbalance wait, ~0.4% to gain |
| KaHIP (team partition, max 37 neighbours, half the inter-node faces) | −0.36% vs scotch `(1 0.5)`, within noise |
| scotch domainWeights / strategies / hierarchical (82 candidates, per-rank cost model) | nothing better than `(4 104)/(1 0.5)` |

## `scripts/` and `wmake-rules/`

PBS job scripts (builds, profiling with perf / VTune / APS / IPM, benchmarks) and the extra
wmake rule files (`OptSPR` and variants) used with the OpenFOAM tree. Paths are specific to
Gadi project storage.

## License

GPL-3.0, as OpenFOAM: `PCGFixed.C` is a modified copy of OpenFOAM's `PCG.C`. See `LICENSE`.
