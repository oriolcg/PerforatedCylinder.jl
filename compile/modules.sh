module load 2022r2 openmpi intel-mkl
export LD_LIBRARY_PATH=$HOME/progs/install/petsc/3.15.4/lib:$LD_LIBRARY_PATH
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/mnt/shared/apps/2022r2/compute/linux-rhel8-skylake_avx512/gcc-8.5.0/openmpi-4.1.1-urzuzcvzrdedifi3mm527t4wgiisuvld
export JULIA_MPIEXEC=srun
export JULIA_PETSC_LIBRARY=$HOME/progs/install/petsc/3.15.4/lib/libpetsc.so
export PerforatedCylinder_MESHES=/scratch/ocolomesgene/tests/julia/PerforatedCylinder.jl/data/meshes
