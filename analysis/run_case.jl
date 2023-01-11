using MPI
MPI.Init()
comm = MPI.COMM_WORLD
using PerforatedCylinder

# Paths
const project_root = ".."
const data_dir = project_root * "/data"
const OVERWRITE = true

# Define cases
nbeta = 3
nalpha = 3
perf_cases = [3,9,27]
porosities = 0.3:(0.7-0.3)/(nbeta-1):0.7
alphas = 0.0:15.0/(nalpha-1):15.0
cases = []
for num_perforations in perf_cases
  for β in porosities
    for α in alphas
      β2 = round(β;digits=2)
      α2 = round(α,digits=2)
      push!(cases,"$num_perforations-$β2-$α2")
    end
  end
end

# Set filenames
np = MPI.Comm_size(comm)
case_id = parse(Int,ENV["CASE_ID"])
testname = cases[case_id]
mesh_file = testname * ".msh"
force_file = testname * ".csv"
output_path = "results_"*testname
if MPI.Comm_rank(comm)==0
  isdir(output_path) || mkdir(output_path)
end
MPI.Barrier(comm)
if MPI.Comm_rank(comm)==0
  println("Testname: $testname")
  println("mesh_file: ",mesh_file)
  println("force_file: ",force_file)
  println("output_path: ",output_path)
  println("Running test case " * testname)
end

# Run case
PerforatedCylinder.main_parallel(np;
  mesh_file=mesh_file,
  force_file=force_file,
  output_path=output_path,
  Δt=0.05,
  tf=100,
  Δtout=0.5,
)

MPI.Finalize()
