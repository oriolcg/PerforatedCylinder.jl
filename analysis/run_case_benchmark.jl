module Run_Case_Parallel
using MPI
MPI.Init()
comm = MPI.COMM_WORLD
using PerforatedCylinder

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const OVERWRITE = true

# Define cases
perf_cases = 3:30
porosities = 0.3:0.02:0.7
# lengths = 5:15
alphas = [0.0]#:15.0/(nalpha-1):15.0
cases = []#["tmp_coarse"]

for num_perforations in perf_cases
  for β in porosities
    α = 0.0
    β2 = round(β;digits=2)
    α2 = round(α,digits=2)
    push!(cases,"$β2-$num_perforations")
  end
end

# Set filenames
np = MPI.Comm_size(comm)
case_id = parse(Int,ENV["CASE_ID"])
testname = cases[case_id]
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
