#!/bin/sh
#
## #SBATCH --account=research-ceg-he
#SBATCH --job-name="PerforatedCylinder"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

source ../compile/modules.sh

mpiexecjl --project=../ -n 1 $HOME/progs/install/julia/1.7.2/bin/julia -J ../PerforatedCylinder.so -O3 --check-bounds=no -e 'using PerforatedCylinder; PerforatedCylinder.main_serial(;mesh_file="test_conformal_mesh.msh",Δt=0.05,tf=0.1,write_vtk=false)'
