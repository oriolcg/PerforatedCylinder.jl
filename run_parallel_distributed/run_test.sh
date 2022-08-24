#!/bin/sh
#
## #SBATCH --account=research-ceg-he
#SBATCH --job-name="CNN_NS"
#SBATCH --partition=compute
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-%j-%4t.out
## #SBATCH --ntasks=480

source ../compile/modules.sh

mpiexecjl --project=../ -n 8 $HOME/progs/install/julia/1.7.2/bin/julia -O3 --check-bounds=no -e 'include("../src/CNN_NS.jl"); CNN_NS.main_parallel_sequential(8; mesh_file="../data/test_conformal_mesh.msh",force_file="forces_tmp.csv",output_path=".")'
