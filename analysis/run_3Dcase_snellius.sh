#!/bin/sh
#
#SBATCH --job-name="perf_cylinder_3D"
#SBATCH --partition=thin
#SBATCH --time=2-00:00:00
#SBATCH -n 32
#SBATCH -o stdout/slurm-%j-%4t.out
#SBATCH -e stdout/slurm-%j-%4t.err

source ../compile/modules_snellius.sh
mpiexecjl --project=../ -n 32 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_3Dcase.jl")'
