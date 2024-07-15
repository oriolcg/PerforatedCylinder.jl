#!/bin/sh
#
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=genoa
#SBATCH --time=1-00:00:00
#SBATCH -n 1
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
echo "case: new_compile"
export CASE_ID="new_compile"
# mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
srun -N1 -n1 -c1 --mem=0 --exact julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no -e 'include("run_full_cylinder_case_serial.jl")'