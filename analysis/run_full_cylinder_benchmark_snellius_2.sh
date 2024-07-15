#!/bin/sh
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=genoa
#SBATCH --time=02:00:00
#SBATCH -n 2
#SBATCH -o stdout-benchmark/slurm-%j-%4t-%n.out
#SBATCH -e stdout-benchmark/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
julia --project=../ -e 'using Pkg; Pkg.instantiate(); Pkg.build("MPI")'
export CASE_ID="new_compile"
echo "Starting case: $CASE_ID"
# mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_benchmark.jl")'
# julia --project=../ -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_benchmark_serial.jl")'
srun -N1 -n2 --mem=0 --exact julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no -e 'include("run_full_cylinder_case_benchmark.jl")'