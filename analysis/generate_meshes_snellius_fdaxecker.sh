#!/bin/bash

#SBATCH --job-name="generate_meshes"
#SBATCH -p genoa
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -o stdout/generate_meshes.out
#SBATCH -e stdout/generate_meshes.err

source ../compile/modules_snellius_fdaxecker.sh
# mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no --color=yes -e 'include("generate_meshes.jl")'
# julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no --color=yes -e 'include("generate_meshes.jl")'
julia --project=../ -O3 --check-bounds=no --color=yes -e 'include("generate_meshes_length.jl")'