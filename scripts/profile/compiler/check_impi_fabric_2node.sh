#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=208,mem=100gb
#PBS -l walltime=00:10:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N impi_fabric_chk
#PBS -j oe

# Confirm what the team's Intel MPI settings actually use between normalsr nodes:
#  1) IB hardware on the compute nodes
#  2) libfabric provider / UCX transports chosen (I_MPI_DEBUG, UCX_PROTO_INFO)
#  3) inter-node PingPong latency and bandwidth vs the IB line rate
#  4) small-message Allreduce on 208 ranks: default vs UCX DC / RC transports
# Same Intel MPI env as run/2node/run.sh (shm:ofi + mlx, scatter pinning off here: 1 rank/node tests).

set -uo pipefail
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/mpi-fabric/check_${PBS_JOBID%%.*}
mkdir -p "$ROOT"; cd "$ROOT"
exec > >(tee job.log) 2>&1
module purge
module load intel-mpi/2021.17.2
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
sort -u "$PBS_NODEFILE" > hosts
echo "== IB on compute nodes"
for h in $(cat hosts); do
    ssh -o BatchMode=yes "$h" 'echo "$(hostname -s): $(for d in /sys/class/infiniband/*; do echo -n "${d##*/} $(cat $d/ports/1/rate) $(cat $d/ports/1/state); "; done)"'
done
echo "== provider/transport selection (I_MPI_DEBUG=5, UCX protocol info)"
I_MPI_DEBUG=5 UCX_PROTO_INFO=y timeout 120 mpirun -np 2 -ppn 1 IMB-MPI1 PingPong -msglog 0:2 > selection.log 2>&1
grep -iE 'libfabric (version|provider)|provider:|of_provider|ucx|transport|mlx5|dc_mlx|rc_mlx|ud_mlx' selection.log | sort -u | head -30
echo "== inter-node PingPong (1 rank/node)"
timeout 180 mpirun -np 2 -ppn 1 IMB-MPI1 PingPong -msglog 0:24 -iter_policy off > pingpong.log 2>&1
awk '/#bytes/{f=1;next} f && NF==4 && $1+0>=0 {printf "  %10d B  lat %8.2f us  bw %9.1f MB/s\n",$1,$3,$4} /^$/{if(f)exit}' pingpong.log | sed -n '1p;3p;11p;17p;21p;25p'
echo "== Allreduce, 208 ranks (104/node), 8 B - 64 KB"
for tls in default dc rc; do
    if [[ $tls == default ]]; then unset UCX_TLS; else export UCX_TLS=self,sm,${tls}_mlx5; fi
    timeout 180 mpirun -np 208 -ppn 104 IMB-MPI1 Allreduce -msglog 3:16 -npmin 208 > allreduce_$tls.log 2>&1
    echo "-- UCX_TLS=${UCX_TLS:-default}"
    awk '/#bytes/{f=1;next} f && NF>=5 && $1+0>0 {printf "  %7d B  avg %8.2f us\n",$1,$5}' allreduce_$tls.log | sed -n '1p;4p;8p;11p;14p'
done
unset UCX_TLS
echo "FABRIC_CHECK_COMPLETE root=$ROOT"
