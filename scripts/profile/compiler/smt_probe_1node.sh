#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=200gb
#PBS -l walltime=00:15:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N smt_probe
#PBS -j oe

# Hyperthreading (SMT) feasibility on one normalsr node:
#   1. does the job's cpuset include the hyperthread siblings (logical CPUs beyond 104)?
#   2. can Intel MPI pin 208 ranks to 208 distinct logical CPUs?
#   3. does a second rank per core add throughput? Same total work over 104 ranks (one per
#      core, team scatter order) and over 208 ranks (both hyperthreads of every core):
#        triad  - streaming, bandwidth-bound (field operations)
#        gather - local-random indirect reads, latency-bound (GaussSeidel on a renumbered mesh)
set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
P=/scratch/ya41/kn3896/openfoam/tools/smtprobe/smtprobe
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/smt/probe_${PBS_JOBID%%.*}
mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
cd "$ROOT"

echo "== allocation"
grep -E 'Cpus_allowed_list' /proc/self/status
echo "nproc=$(nproc)"
lscpu | grep -E '^(CPU\(s\)|Thread\(s\) per core|Core\(s\) per socket|Socket\(s\)|NUMA node\(s\))'
cat /sys/fs/cgroup/cpuset.cpus.effective 2>/dev/null || cat /sys/fs/cgroup/cpuset/cpuset.cpus 2>/dev/null || find /sys/fs/cgroup -maxdepth 4 -name 'cpuset.cpus*' 2>/dev/null | head -3

set +eu
module purge
module load intel-mpi/2021.17.2
set -eu
export I_MPI_FABRICS=shm I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0 I_MPI_PIN=1 I_MPI_HYDRA_TOPOLIB=hwloc
P104=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$TEAM/run/1node/intel-scatter-cpus.awk")
P208=$(for c in ${P104//,/ }; do s=$(cut -d, -f1-2 /sys/devices/system/cpu/cpu$c/topology/thread_siblings_list | tr '-' ','); echo "$s"; done | tr ',' '\n' | awk '!seen[$0]++' | paste -sd,)
echo "pin104: $(tr ',' '\n' <<< "$P104" | wc -l) cpus; pin208: $(tr ',' '\n' <<< "$P208" | wc -l) cpus"
echo "thread siblings of cpu ${P104%%,*}: $(cat /sys/devices/system/cpu/cpu${P104%%,*}/topology/thread_siblings_list)"

echo "== placement check (208 ranks)"
I_MPI_PIN_PROCESSOR_LIST=$P208 mpirun -np 208 -ppn 208 sh -c 'taskset -cp $$ | awk "{print \$NF}"' > pin208.txt 2>&1 || true
echo "ranks reporting: $(wc -l < pin208.txt), distinct single-cpu masks: $(grep -cE '^[0-9]+$' pin208.txt), unique cpus: $(grep -E '^[0-9]+$' pin208.txt | sort -u | wc -l), max cpu: $(grep -E '^[0-9]+$' pin208.txt | sort -n | tail -1)"

bench() {   # $1 mode, $2 ranks, $3 elements per rank, $4 reps, $5 pin list
    local f=bench_$1_$2.txt
    I_MPI_PIN_PROCESSOR_LIST=$5 mpirun -np $2 -ppn $2 "$P" "$1" "$3" "$4" > "$f" 2>&1
    awk -v m=$1 -v r=$2 '{t=$1; if (t>max) max=t; s+=t} END {printf "%-6s ranks=%3d  makespan %.3f s  mean %.3f s\n", m, r, max, s/NR}' "$f"
}
echo "== throughput: same total work (104 x 4M = 208 x 2M elements)"
for rep in 1 2; do
    bench triad  104 4000000 30 "$P104"
    bench triad  208 2000000 30 "$P208"
    bench gather 104 4000000 30 "$P104"
    bench gather 208 2000000 30 "$P208"
done
echo "SMT_PROBE_COMPLETE root=$ROOT"
