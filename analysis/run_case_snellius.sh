#!/bin/sh
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=thin
#SBATCH --time=2-00:00:00
#SBATCH -n 12
#SBATCH -o stdout/slurm-%j-%4t.out
#SBATCH -e stdout/slurm-%j-%4t.err

source ../compile/modules_snellius.sh
export CASE_ID=$1
echo "Starting case: $CASE_ID"
mpiexecjl --project=../ -n 12 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case.jl")'
