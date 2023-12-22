module load 2022 OpenMPI/4.1.4-GCC-11.3.0 imkl/2022.1.0 Julia/1.8.2-linux-x86_64
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/sw/arch/RHEL8/EB_production/2022/software/OpenMPI/4.1.4-GCC-11.3.0
export JULIA_PETSC_LIBRARY=/gpfs/home3/colomeso/progs/install/petsc/3.18/lib/libpetsc
export PerforatedCylinder_MESHES=/gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/meshes
export PerforatedCylinder_FORCES=/gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/forces 
export PerforatedCylinder_VTKs=/gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/VTKs