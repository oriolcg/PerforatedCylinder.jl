#!/bin/bash

#SBATCH --job-name="cnn-NS"
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G

source modules.sh

mpiexecjl --project=../ -n 1 $HOME/progs/install/julia/1.7.2/bin/julia -O3 --check-bounds=no --color=yes compile.jl

