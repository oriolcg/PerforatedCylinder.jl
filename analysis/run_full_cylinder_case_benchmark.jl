module Run_Case_Parallel
using MPI
MPI.Init()
comm = MPI.COMM_WORLD
using PerforatedCylinder

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const OVERWRITE = true


# Set filenames
np = MPI.Comm_size(comm)
testname = ENV["CASE_ID"]
mesh_file = testname * ".msh"
force_file = "parallel_" * testname * ".csv"
vtks_path = ENV["PerforatedCylinder_VTKs"]
output_path = joinpath(vtks_path,"parallel_results_"*testname*"-$np")
if isdir(output_path)
  if MPI.Comm_rank(comm)==0
    println("Existing case. Exiting execution without computing.")
  end
  MPI.Finalize()
  return nothing
end
MPI.Barrier(comm)
if MPI.Comm_rank(comm)==0
   mkdir(output_path)
   println("Testname: $testname")
   println("mesh_file: ",mesh_file)
   println("force_file: ",force_file)
   println("output_path: ",output_path)
   println("Running test case " * testname)
end
MPI.Barrier(comm)

# Run case
PerforatedCylinder.main_parallel(np;
  mesh_file=mesh_file,
  force_file=force_file,
  output_path=output_path,
  Δt=0.01,
  tf=0.5,
  Δtout=0.05,
)

MPI.Finalize()
end
