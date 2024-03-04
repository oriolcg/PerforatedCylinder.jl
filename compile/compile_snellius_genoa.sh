#!/bin/bash

#SBATCH --job-name="compile_PerforatedCylinder"
#SBATCH -p genoa
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -o stdout_genoa
#SBATCH -e stderr_genoa

source modules_snellius.sh
julia --project=../ -e 'using Pkg; Pkg.build("MPI")'
mpiexecjl --project=../ -n 1 julia -O3 --check-bounds=no --color=yes compile_genoa.jl

