#!/bin/bash

#SBATCH --job-name="generate_meshes"
#SBATCH -p thin
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -o stdout/generate_meshes.out
#SBATCH -e stdout/generate_meshes.err

source ../compile/modules_snellius.sh
mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no --color=yes -e 'include("generate_meshes.jl")'