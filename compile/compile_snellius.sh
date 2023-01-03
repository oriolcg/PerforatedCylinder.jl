#!/bin/bash

#SBATCH --job-name="compile_CNN_NS"
#SBATCH -p thin
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -o stdout
#SBATCH -e stderr

source modules_snellius.sh

mpiexecjl --project=../ -n 1 julia -O3 --check-bounds=no --color=yes compile.jl

