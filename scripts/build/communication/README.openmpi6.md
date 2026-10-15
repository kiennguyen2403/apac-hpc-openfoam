# Open MPI 6 experiment

Original source: `/home/565/kn3896/src/openmpi-6.0.0rc2`, cloned recursively from
https://github.com/open-mpi/ompi.git at tag `v6.0.0rc2`, commit
`733f33ec094d58b73d38c1c948265cf4aa211e3a`. This is a release candidate,
not a final release. The source includes tuned Allreduce algorithm 7,
`allgather_reduce`.

Build job `180341537.gadi-pbs` completed successfully (12m32s, 13.37 SU).
The clean source checkout, including its Git metadata and bundled dependency
sources, was removed at the user's request after installation. The installed
binaries, libraries, headers, build script and build logs were retained.

To restore the source for a future rebuild:

```bash
git clone --branch v6.0.0rc2 --depth 1 --recurse-submodules --shallow-submodules \
    https://github.com/open-mpi/ompi.git /home/565/kn3896/src/openmpi-6.0.0rc2
```

The source is not needed to run the installed MPI. The saved build script is
`build_openmpi6.sh`; it deliberately refuses to overwrite the existing install.

The job requests 32 CPUs, 64 GB RAM, 12 GB job-local storage and one hour
on Gadi's normalsr (Sapphire Rapids) queue. Compilation uses PBS job-local
storage to avoid scratch file-quota pressure. The installed
files still count towards scratch quota.

It builds with ICX 2025.3.2 and
`-Ofast -march=sapphirerapids -mprefer-vector-width=512 -fPIC` for both C and C++,
without debug instrumentation, Fortran bindings or OSHMEM. These flags
target Sapphire Rapids and relax floating-point rules; validate numerical
results before relying on this build. ICX accepts `-Ofast` but emits a
deprecation warning. The resulting binaries require compatible CPUs.
It reuses site UCX 1.18.1, hwloc 2.12.1,
libevent 2.1.12, PMIx 6.0.0 and PRRTE 4.0.0. The site's PRRTE contains
the PBS/TM launcher. Autoconf 2.72 and Python 3.11.7 are build tools.

Install: `/scratch/ya41/kn3896/openfoam/ThirdParty-v2512/platforms/linux64Icx/openmpi-6.0.0rc2`

Logs: `/scratch/ya41/kn3896/openfoam/logs/build/openmpi6_<PBS_JOBID>/`

The script stops if that install directory already exists. On failure,
inspect the stage logs before retrying; a failed installation may leave
a partial install directory. Build intermediates disappear with job-local
storage when the job ends.

Post-build checks verify the algorithm listing, UCX PML and PBS launcher
component. They do not launch MPI ranks or prove multi-node correctness.
Existing MPI installations, OpenFOAM builds and submission scripts are
unchanged. Rebuilding an isolated OpenFOAM MPI variant and running tiny
multi-node tests are separate next steps. Algorithm 7 is not forced by
this build script; its availability does not guarantee automatic selection
or a performance gain on Gadi.
