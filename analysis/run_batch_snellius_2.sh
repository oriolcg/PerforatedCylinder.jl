#!/bin/sh
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=genoa
#SBATCH --time=02:00:00
#SBATCH -n 2
#SBATCH -o stdout-benchmark/slurm-%j-%4t-%n.out
#SBATCH -e stdout-benchmark/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
INITIAL_CASE=1
FINAL_CASE=1
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_benchmark.jl")'
    # julia --project=../ -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_benchmark_serial.jl")'
    srun -N1 -n2 --mem=0 --exact julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no -e 'include("run_case_benchmark.jl")' &
done
wait