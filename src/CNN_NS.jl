module CNN_NS

using Gridap
using Gridap.FESpaces: zero_free_values, interpolate!
using Gridap.Fields: meas
using GridapGmsh
using GridapDistributed
using GridapDistributed: DistributedTriangulation, DistributedCellField
using GridapPETSc
using GridapPETSc: PETSC
using PartitionedArrays
using SparseMatricesCSR

using CSV
using DataFrames

include("NavierStokesSerial.jl")
include("NavierStokesParallel.jl")

function main_parallel(np;
  mesh_file="test_conformal_mesh_coarse.msh",
  force_file="forces.csv",
  output_path="tmp",
  Δt=0.5,
  tf=1.0,
  Δtout=0.5)
  prun(mpi,np) do parts
    options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-12 -snes_atol 0.0 -snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0"
    GridapPETSc.with(args=split(options)) do
      run_test_parallel(parts,mesh_file,force_file,output_path,Δt,tf,Δtout)
    end
  end
end

function main_parallel_sequential(np;
  mesh_file="test_conformal_mesh_coarse.msh",
  force_file="forces.csv",
  output_path="tmp",
  Δt=0.5,
  tf=1.0,
  Δtout=0.5)
  prun(sequential,np) do parts
    options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-12 -snes_atol 0.0 -snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0"
    GridapPETSc.with(args=split(options)) do
      run_test_parallel(parts,mesh_file,force_file,output_path,Δt,tf,Δtout)
    end
  end
end
# function main_parallel(np;mesh_file="test_conformal_mesh_coarse.msh",Δt=0.5,tf=1.0)
#   prun(mpi,np) do parts
#     options = "-snes_type newtonls -snes_linesearch_type basic  -snes_linesearch_damping 1.0 -snes_rtol 1.0e-12 -snes_atol 0.0 -snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
#     GridapPETSc.with(args=split(options)) do
#       FD, FL = run_test_parallel(parts,mesh_file,Δt,tf)

#       println(FD)
#       println(FL)

#       forces_path = ENV["CNN_NS_FORCES"]
#       fname = splitext(mesh_file)[1]
#       forces_file = joinpath(forces_path,"$fname.csv")
#       T = t₀:Δt:tf
#       CSV.write(forces_file, DataFrame(t = T, FD = FD, FL = FL))
#     end
#   end
# end

function main_serial(;
  mesh_file="test_conformal_mesh_coarse.msh",
  force_file="forces.csv",
  output_path="tmp",
  Δt=0.5,
  tf=1.0,
  Δtout=0.5)
  run_test_serial(mesh_file,force_file,output_path,Δt,tf,Δtout)

  # forces_path = ENV["CNN_NS_FORCES"]
  # fname = splitext(mesh_file)[1]
  # forces_file = joinpath(forces_path,"$fname.csv")
  # T = t₀:Δt:tf
  # CSV.write(forces_file, DataFrame(t = T, FD = FD, FL = FL))
end

end
