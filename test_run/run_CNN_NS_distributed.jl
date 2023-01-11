using MPI
MPI.Init()
comm = MPI.COMM_WORLD

using PerforatedCylinder
using CSV
using DataFrames

# Paths
const project_root = ".."
const data_dir = project_root * "/data"
const forces_dir = data_dir * "/forces/"
const meshes_dir = data_dir * "/meshes/"
const OVERWRITE = true

testname = "ID_" * string(MPI.Comm_rank(comm)+1, pad=6)
mesh_file = meshes_dir * testname * ".msh"
force_file = forces_dir * testname * ".csv"
output_path = joinpath(pwd(),"results_"*testname)
println("Running test case " * testname)
isdir(output_path) || mkdir(output_path)
FD,FL = PerforatedCylinder.main_serial(;
  mesh_file=mesh_file,
  Δt=0.05,
  tf=200.0,
  write_vtk=true,
  output_path=output_path,
  Δtout=0.5,
)
CSV.write(force_file, DataFrame(FD = FD, FL = FL))

#=

isdir(forces_dir) || mkdir(forces_dir)

function main()
  # %% Processing
  println("Start main function")

  for (root, dirs, files) in walkdir(meshes_dir)

      println("Start parallel work")
      println(root," ",dirs," ",files)

      @assert(length(files) >= nworkers())
      @sync @distributed for file in files

        fname = splitext(file)[1]
        id = parse(Int64, fname[4:end])
        forces_file = forces_dir * fname * ".csv"
        println(forces_file)

        """Skip work if we already have the files and don't want to overwrite them.

        This will save much time when restarting a previously interrupted run.
        Under the condition that previously generated samples are all of good quality.

        TODO this causes problems with incomplete csv results. Need to fix.
        As csv will not be written when interrupted.
        Perhaps, also check if the cfd results exists in the dataframe before skipping.
        """
        if OVERWRITE || sample_doesnt_already_exists([forces_file])
            mesh_path = joinpath(root, file)
            testname = replace(file,".msh" =>"")
            output_path = joinpath(pwd(),"results_"*testname)
            println("Running test case " * testname)
            isdir(output_path) || mkdir(output_path)
            FD,FL = PerforatedCylinder.main_serial(;mesh_file=mesh_path,Δt=0.05,tf=200,write_vtk=false,output_path)
            CSV.write(forces_file, DataFrame(FD = FD, FL = FL))
        end
      end
  end
  println("Finished all files")
end

main()

rmprocs(workers())
=#
MPI.Finalize()
