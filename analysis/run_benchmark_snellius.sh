#!/bin/sh
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=thin
#SBATCH --time=2-00:00:00
#SBATCH -n ${NP}
#SBATCH -o stdout-benchmark/slurm-%j-%4t-%n.out
#SBATCH -e stdout-benchmark/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
export CASE_ID=1
echo "Starting case: $CASE_ID"
mpiexecjl --project=../ -n $1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_benchmark.jl")'
