using MPI
MPI.Init()
comm = MPI.COMM_WORLD
initial_case = parse(Int,ENV["INITIAL_CASE"])
local_case = parse(Int,ENV["LOCAL_CASE"])

using CNN_NS

# Paths
const project_root = ".."
const data_dir = project_root * "/data"
const forces_dir = data_dir * "/forces/"
const meshes_dir = data_dir * "/meshes/"
const OVERWRITE = true

testname = "ID_" * string((((initial_case-1)+(local_case-1))*48 + (MPI.Comm_rank(comm)+1)), pad=6)
println("Testname: $testname")
mesh_file = meshes_dir * testname * ".msh"
force_file = forces_dir * testname * ".csv"
output_path = joinpath(pwd(),"results_"*testname)
println("Running test case " * testname)
isdir(output_path) || mkdir(output_path)
stdout_file = joinpath(output_path,"stdout")
stderr_file = joinpath(output_path,"stderr")
redirect_stdio(;stdout=stdout,stderr=stderr) do
  CNN_NS.main_serial(;
    mesh_file=mesh_file,
    force_file=force_file,
    output_path=output_path,
    Δt=0.1,
    tf=400.0,
    Δtout=0.5,
  )
end
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
            FD,FL = CNN_NS.main_serial(;mesh_file=mesh_path,Δt=0.05,tf=200,write_vtk=false,output_path)
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
