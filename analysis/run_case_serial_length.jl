module Run_Case_Serial
using PerforatedCylinder
using TimerOutputs

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
# cases = [ "0.42-13",
#           "0.42-19",
#           "0.44-9",
#           "0.46-15",
#           "0.5-28",
#           "0.68-22",
#           "0.68-30",
#           "0.68-6",
#           "0.33-11",
#           "0.33-12",
#           "0.33-13",
#           "0.33-14",
#           "0.33-15",
#           "0.33-16",
#           "0.33-17",
#           "0.33-18",
#           "0.33-19",
#           "0.33-20",
#           "0.33-21",
#           "0.33-22",
#           "0.33-23",
#           "0.33-24",
#           "0.33-25",
#           "0.33-26",
#           "0.33-27",
#           "0.33-28",
#           "0.33-29",
#           "0.33-30",
#           "0.33-31",
#           "0.33-32",
#           "0.33-33",
#           "0.33-34",
#           "0.33-35",
#           "0.33-36",
#           "0.33-37",
#           "0.33-38",
#           "0.33-39",
#           "0.33-40",
#           "0.33-41",
#           "0.33-42",
#           "0.33-43",
#           "0.33-44",
#           "0.33-45",
#           "0.33-46",
#           "0.33-47",
#           "0.33-48",
#           "0.33-49",
#           "0.33-50",
#           "0.33-51",
#           "0.33-52",
#           "0.33-53",
#           "0.33-54",
#           "0.33-55",
#           "0.33-56",
#           "0.33-57"]
# for i in 1:3
#   for num_perforations in perf_cases
#       for β in porosities
#           β2 = round(β, digits=2)
#           push!(cases, "$β2-$num_perforations")
#       end
#   end
# end

# const to = TimerOutput()

# @timeit to "coarse" begin

# Set filenames
# for testname in cases[1:10]
# time_steps = [0.2, 0.1, 0.05]
# for Δt in time_steps
  testname = cases[parse(Int,ENV["CASE_ID"])]
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
