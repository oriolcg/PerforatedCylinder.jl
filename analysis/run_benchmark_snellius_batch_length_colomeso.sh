#!/bin/sh
#
#SBATCH --mem=140G
#SBATCH --job-name="perf_cylinder"
#SBATCH --partition=genoa
#SBATCH --time=1-00:00:00
#SBATCH -n 80
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err

# source ../compile/modules_snellius.sh
source ../compile/modules_snellius.sh
INITIAL_CASE=1
FINAL_CASE=16
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
    srun -N1 -n1 -c5 --mem=0 --exact julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no -e 'include("run_case_serial_length.jl")' &
done
wait
