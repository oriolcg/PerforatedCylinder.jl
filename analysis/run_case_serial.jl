module Run_Case_Serial
using PerforatedCylinder

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const OVERWRITE = true
ENV["PerforatedCylinder_VTKs"] = data_dir * "/VTKs/"
ENV["PerforatedCylinder_FORCES"] = data_dir * "/forces/"
ENV["PerforatedCylinder_MESHES"] = data_dir * "/meshes/"

# Define cases
nbeta = 20
nalpha = 1
nperfs = 40
perf_cases = [3]#[3,9,27]
porosities = [0.3]#0.3:(0.7-0.3)/(nbeta-1):0.7
alphas = [0.0]#:15.0/(nalpha-1):15.0
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
case_id = 1
testname = cases[case_id]
mesh_file = testname * ".msh"
force_file = testname * ".csv"
vtks_path = ENV["PerforatedCylinder_VTKs"]
output_path = joinpath(vtks_path,"results_"*testname)
if isdir(output_path)
  println("Existing case. Exiting execution without computing.")
  return nothing
else
  mkdir(output_path)
end
println("Testname: $testname")
println("mesh_file: ",mesh_file)
println("force_file: ",force_file)
println("output_path: ",output_path)
println("Running test case " * testname)

# Run case
PerforatedCylinder.main_serial(mesh_file,
  force_file,output_path,0.05,1.0,0.0
)
end
