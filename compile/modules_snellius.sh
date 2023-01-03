module load 2021 OpenMPI/4.1.1-GCC-10.3.0 imkl/2021.4.0 Julia/1.7.3-linux-x86_64
export LD_LIBRARY_PATH=$HOME/progs/install/petsc/3.15.4/lib:$LD_LIBRARY_PATH
export PATH=$HOME/.julia/bin:$PATH
export JULIA_DEPOT_PATH=$HOME/.julia
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/sw/arch/Centos8/EB_production/2021/software/OpenMPI/4.1.1-GCC-10.3.0
export JULIA_MPIEXEC=srun
export JULIA_PETSC_LIBRARY=$HOME/progs/install/petsc/3.15.4/lib/libpetsc.so
export CNN_NS_MESHES=/gpfs/scratch1/nodespecific/int1/colomeso/tests/CNN_NS.jl/data/meshes
export CNN_NS_RESULTS=/gpfs/scratch1/nodespecific/int1/colomeso/tests/CNN_NS.jl 