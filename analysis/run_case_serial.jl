module Run_Case_Serial
using PerforatedCylinder

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const OVERWRITE = true

# Define cases
nbeta = 21
nalpha = 1
nperfs = 40
perf_cases = 3:30
porosities = 0.3:0.02:0.7
alphas = [0.0]#:15.0/(nalpha-1):15.0
cases = []#["tmp_coarse"]
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
# for testname in cases[1:10]
  testname = cases[parse(Int,ENV["CASE_ID"])]
  mesh_file = testname * ".msh"
  force_file = testname * ".csv"
  vtks_path = ENV["PerforatedCylinder_VTKs"]
  output_path = joinpath(vtks_path,"results_"*testname)
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
    Δt = 0.05
    tf = 1.0
    Δtout = 0.0
  end
  PerforatedCylinder.main_serial(mesh_file=mesh_file,
    force_file=force_file,
    output_path=output_path,
    Δt=Δt,
    tf=tf,
    Δtout=Δtout)
# end
end
