# OpenFOAM v2512 performance work — APAC HPC-AI 2026 (occDrivAerStaticMesh 65M)

Tools, solver code and job scripts from tuning `simpleFoam` on the 65M-cell
occDrivAerStaticMesh case on NCI Gadi (Sapphire Rapids, `normalsr`), OpenFOAM v2512
built with Intel oneAPI ICX 2025.3 and Intel MPI 2021.17.

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

Enable in a case:

```
// system/controlDict
libs ("<path>/libfixedIterSolvers.so");

// system/fvSolution: smoothSolver -> smoothSolverFixed, coarsest-level PCG -> PCGFixed
```

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

## `scripts/` and `wmake-rules/`

PBS job scripts (builds, profiling with perf / VTune / APS / IPM, benchmarks) and the extra
wmake rule files (`OptSPR` and variants) used with the OpenFOAM tree. Paths are specific to
Gadi project storage.

## License

GPL-3.0, as OpenFOAM: `PCGFixed.C` is a modified copy of OpenFOAM's `PCG.C`. See `LICENSE`.
