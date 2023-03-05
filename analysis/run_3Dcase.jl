using MPI
MPI.Init()
comm = MPI.COMM_WORLD
using PerforatedCylinder

# Paths
const project_root = ".."
const data_dir = project_root * "/data"
const OVERWRITE = true

# Set filenames
np = MPI.Comm_size(comm)
testname = "3D_monopile_coarse"
mesh_file = "3D_monopile_coarse.msh"
force_file = "3D_monopile_coarse.csv"
vtks_path = ENV["PerforatedCylinder_VTKs"]
output_path = joinpath(vtks_path,"results_"*testname)
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
  Δt=0.05,
  tf=1.0,
  Δtout=0.5,
)

MPI.Finalize()
