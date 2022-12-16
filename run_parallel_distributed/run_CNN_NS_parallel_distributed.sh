#!/bin/sh
#
## #SBATCH --account=research-ceg-he
#SBATCH --job-name="CNN_NS"
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
## #SBATCH --exclusive
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-%j-%4t.out
## #SBATCH --ntasks=480

source ../compile/modules.sh

mpiexecjl --project=../ -n 24 $HOME/progs/install/julia/1.7.2/bin/julia -J ../CNN_NS_parallel.so -O3 --check-bounds=no -e 'include("run_CNN_NS_parallel_distributed.jl")'
