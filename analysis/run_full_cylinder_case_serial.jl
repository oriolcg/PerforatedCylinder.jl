module Run_Full_Cylinder_Case_Serial
using PerforatedCylinder
using TimerOutputs

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const OVERWRITE = true
testname = ENV["CASE_ID"]#cases[parse(Int,ENV["CASE_ID"])]
mesh_file = testname * ".msh"
force_file = testname * ".csv"
vtks_path = ENV["PerforatedCylinder_VTKs"]
output_path = joinpath(vtks_path,"serial_results_"*testname)
# output_path = joinpath(vtks_path,"results_"*testname)
if isdir(output_path)
  println("Existing case. Exiting execution without computing.")
else
  mkdir(output_path)
end
println("Testname: $testname")
println("mesh_file: ",mesh_file)
println("force_file: ",force_file)
println("output_path: ",output_path)
println("Running test case " * testname)

# Run case
if testname == "tmp_coarse"
  Δt = 0.05
  tf = Δt
  Δtout = 0.0
else
  Δt = 0.01
  # Δt = Δt
  tf = 50.0
  Δtout = 0.05
end
PerforatedCylinder.main_serial(mesh_file=mesh_file,
  force_file=force_file,
  output_path=output_path,
  Δt=Δt,
  tf=tf,
  Δtout=Δtout)
# end
# end

# show(to)
end
