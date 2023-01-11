module load 2022 OpenMPI/4.1.4-GCC-11.3.0 imkl/2022.1.0 Julia/1.7.3-linux-x86_64
export LD_LIBRARY_PATH=$HOME/progs/install/petsc/3.15.4/lib:$LD_LIBRARY_PATH
export PATH=$HOME/.julia/bin:$PATH
export JULIA_DEPOT_PATH=$HOME/.julia
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/sw/arch/Centos8/EB_production/2022/software/OpenMPI/4.1.4-GCC-11.3.0
export JULIA_MPIEXEC=srun
export JULIA_PETSC_LIBRARY=$HOME/progs/install/petsc/3.15.4/lib/libpetsc.so
export PerforatedCylinder_MESHES=/gpfs/scratch1/nodespecific/int1/colomeso/tests/PerforatedCylinder.jl/data/meshes
export PerforatedCylinder_FORCES=/gpfs/scratch1/nodespecific/int1/colomeso/tests/PerforatedCylinder.jl/data/forces 
export PerforatedCylinder_RESULTS=/gpfs/scratch1/nodespecific/int1/colomeso/tests/PerforatedCylinder.jl 