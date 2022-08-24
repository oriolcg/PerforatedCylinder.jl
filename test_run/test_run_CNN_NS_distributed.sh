#!/bin/sh
#
## #SBATCH --account=research-ceg-he
#SBATCH --job-name="CNN_NS"
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=slurm-%j-%4t.out

source ../compile/modules.sh

mpiexecjl --project=../ -n 1 $HOME/progs/install/julia/1.7.2/bin/julia -J ../CNN_NS.so -O3 --check-bounds=no -e 'include("run_CNN_NS_distributed.jl")'
