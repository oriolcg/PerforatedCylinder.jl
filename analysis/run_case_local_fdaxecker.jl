module Local_run
using PerforatedCylinder
using TimerOutputs
using ProfileView

# Paths
const project_root = joinpath(@__DIR__,"..")
const data_dir = project_root * "/data"
const mesh_dir = data_dir * "/meshes"
const forces_dir = data_dir * "/forces"
const output_dir = data_dir * "/VTKs"
const OVERWRITE = true
ENV["PerforatedCylinder_FORCES"] = forces_dir
ENV["PerforatedCylinder_MESHES"] = mesh_dir
ENV["PerforatedCylinder_VTKs"] = output_dir

const to = TimerOutput()

println("COARSE")
@timeit to "coarse" begin
# Run case
PerforatedCylinder.main_serial(
  mesh_file="tmp_coarse.msh",
  force_file="tmp_coarse.csv",
  output_path=output_dir,
  Δt=0.05,
  tf=0.05,
  Δtout=0.05
)
end
println("COARSE 2")
@timeit to "coarse2" begin
  # Run case
  @profview PerforatedCylinder.main_serial(
    mesh_file="tmp_coarse.msh",
    force_file="tmp_coarse.csv",
    output_path=output_dir,
    Δt=0.05,
    tf=0.05,
    Δtout=0.05
  )
end
# println("COARSE BP")
# @timeit to "coarse_bp" begin
#   # Run case
#   PerforatedCylinder.main_serial_bp(
#     mesh_file="tmp_coarse.msh",
#     force_file="tmp_coarse.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.05,
#     Δtout=0.05
#   )
#   end

# println("COARSE BP 2")
#   @timeit to "coarse_bp2" begin
#     # Run case
#     PerforatedCylinder.main_serial_bp(
#       mesh_file="tmp_coarse.msh",
#       force_file="tmp_coarse.csv",
#       output_path=output_dir,
#       Δt=0.05,
#       tf=0.05,
#       Δtout=0.05
#     )
#   end
# @timeit to "coarse_long" begin
#   # Run case
#   PerforatedCylinder.main_serial(
#     mesh_file="tmp_coarse.msh",
#     force_file="tmp_coarse.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.5,
#     Δtout=0.05
#   )
#   end

# @timeit to "fine" begin
#   # Run case
#   PerforatedCylinder.main_serial(
#     mesh_file="H=2.0-3-0.3-30.0.msh",
#     force_file="H=2.0-3-0.3-30.0.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.01,
#     Δtout=0.05
#   )
# end
# @timeit to "fine_bp" begin
#   # Run case
#   PerforatedCylinder.main_serial_bp(
#     mesh_file="H=2.0-3-0.3-30.0.msh",
#     force_file="H=2.0-3-0.3-30.0.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.01,
#     Δtout=0.05
#   )
# end
# @timeit to "fine long" begin
#   # Run case
#   PerforatedCylinder.main_serial(
#     mesh_file="H=2.0-3-0.3-30.0.msh",
#     force_file="H=2.0-3-0.3-30.0.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.5,
#     Δtout=0.05
#   )
# end
# @timeit to "fine long bp" begin
#   # Run case
#   PerforatedCylinder.main_serial_bp(
#     mesh_file="H=2.0-3-0.3-30.0.msh",
#     force_file="H=2.0-3-0.3-30.0.csv",
#     output_path=output_dir,
#     Δt=0.05,
#     tf=0.5,
#     Δtout=0.05
#   )
# end

show(to)
# for mesh_file_i in mesh_dir:
#   PerforatedCylinder.main_serial(
#   mesh_file=mesh_file_i,
#   force_file="$mesh_file_i.csv",
#   output_path=output_dir,
#   Δt=0.05,
#   tf=30,
#   Δtout=0.1
# )

end
