module PerforatedCylinder

using Gridap
using Gridap.FESpaces: zero_free_values, interpolate!
using Gridap.Fields: meas
using Gridap.ODEs: SDIRK_3_3, EXRK_RungeKutta_4_4, SDIRK_Midpoint_1_2
using Gridap.ODEs: DIMRungeKutta, EXRungeKutta
using GridapGmsh: gmsh, GmshDiscreteModel
using GridapDistributed
using GridapDistributed: DistributedTriangulation, DistributedCellField
using GridapPETSc
using GridapPETSc: PETSC
using PartitionedArrays
using SparseMatricesCSR
using LineSearches: Static, BackTracking
# using GridapSolvers
# using GridapSolvers.LinearSolvers, GridapSolvers.NonlinearSolvers
# using GridapSolvers.BlockSolvers: NonlinearSystemBlock, BlockTriangularSolver

using CSV
using DataFrames

# include("NavierStokesSerial.jl")
include("NavierStokesSerial_coupled.jl")
include("NavierStokesParallel.jl")
include("mesh_generation.jl")
include("mesh_generation_length.jl")

function generate_meshes(nperf=1,nβ=1,nα=1)
  # Create cases
  perforations = 3:3+nperf
  porosities = 0.2:(0.8-0.2)/(nβ):0.8
  angles = 0:15/(nα):15
  for num_perforations in 3:3+nperf
    for β in porosities
      for α in angles
        create_mesh(num_perforations = num_perforations, β = β, α = α)
      end
    end
  end

end

options_pcasm = """
-snes_type newtonls
-snes_linesearch_type basic
-snes_linesearch_damping 1.0
-snes_rtol 1.0e-6
-snes_atol 1.0e-9
-snes_monitor
-snes_converged_reason
-ksp_type gmres
-ksp_rtol 1.0e-4
-ksp_atol 1.0e-14
-ksp_gmres_restart 5000
-ksp_monitor
-pc_type asm
-sub_ksp_type preonly
-sub_pc_type lu
-mm_ksp_type cg
-mm_ksp_monitor
-mm_ksp_rtol 1.0e-4
-mm_pc_type jacobi
"""

options_fieldsplit = "-snes_type newtonls \
-snes_linesearch_type basic \
-snes_linesearch_damping 1.0 \
-snes_rtol 1.0e-12 \
-snes_atol 0.0 \
-snes_monitor \
-ksp_error_if_not_converged true \
-ksp_converged_reason \
-ksp_type fgmres \
-ksp_rtol 1.0e-6 \
-ksp_atol 1.0e-8 \
-ksp_monitor \
-ksp_monitor_true_residual \
-pc_use_amat \
-pc_type fieldsplit \
-pc_fieldsplit_block_size 3 \
-pc_fieldsplit_0_fields 0,1 \
-pc_fieldsplit_1_fields 2 \
-pc_fieldsplit_type schur \
-pc_fieldsplit_schur_fact_type full \
-pc_fieldsplit_schur_precondition selfp "
#=
-pc_fieldsplit_off_diag_use_amat \
-fieldsplit_0_ksp_type gmres \
-fieldsplit_0_ksp_rtol 1.0e-06 \
-fieldsplit_0_ksp_atol 0.0 \
-fieldsplit_0_ksp_monitor \
-fieldsplit_0_pc_type gamg \
-fieldsplit_0_pc_gamg_type agg \
-fieldsplit_0_pc_gamg_est_ksp_type gmres \
-fieldsplit_0_pc_gamg_agg_nsmooths 1 \
-fieldsplit_0_mg_coarse_sub_pc_type cholesky \
-fieldsplit_0_mg_coarse_sub_pc_factor_mat_ordering_type nd \
-fieldsplit_1_mat_schur_complement_ainv_type lump \
-fieldsplit_1_ksp_type gmres \
-fieldsplit_1_ksp_rtol 1.0e-6 \
-fieldsplit_1_ksp_atol 0.0 \
-fieldsplit_0_ksp_monitor "
=#

options_mumps = "-snes_type newtonls \
-snes_linesearch_type basic  \
-snes_linesearch_damping 1.0 \
-snes_rtol 1.0e-6 \
-snes_atol 1.0e-8 \
-snes_max_it 20 \
-ksp_error_if_not_converged true \
-ksp_converged_reason -ksp_type preonly \
-pc_type lu \
-pc_factor_mat_solver_type mumps \
-mat_mumps_icntl_7 0 \
-mat_mumps_icntl_14 500000"


function main_parallel(np;
  mesh_file="tmp_coarse.msh",
  force_file="forces.csv",
  output_path="tmp",
  Δt=0.5,
  tf=1.0,
  Δtout=0.5)
  current_path = pwd()
  cd(output_path)
  with_mpi() do distribute
    options = options_mumps
    ranks = distribute_with_mpi(LinearIndices((np,)))
    GridapPETSc.with(args=split(options)) do
      run_test_parallel(ranks,mesh_file,force_file,Δt,tf,Δtout)
    end
  end
  cd(current_path)
end

function main_serial(;mesh_file="tmp_coarse.msh",
  force_file="forces.csv",
  output_path="tmp",
  Δt=0.5,
  tf=1.0,
  Δtout=0.5)
  println("Running serial test: $(mesh_file), $(force_file), $Δt, $tf, $Δtout")
  current_path = pwd()
  cd(output_path)
  run_test_serial(mesh_file,force_file,Δt,tf,Δtout)
  cd(current_path)
end

end
