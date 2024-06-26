#!/bin/sh
#
#SBATCH --account=research-ceg-he
#SBATCH --job-name="PerforatedCylinder"
#SBATCH --partition=compute
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=3G
#SBATCH --output=slurm-%j-%4t.out
## #SBATCH --ntasks=480

source ../compile/modules.sh

mpiexecjl --project=../ -n 48 $HOME/progs/install/julia/1.7.2/bin/julia -J ../PerforatedCylinder_serial.so -O3 --check-bounds=no -e 'include("run_PerforatedCylinder_distributed.jl")'
