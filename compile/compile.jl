using PackageCompiler
create_sysimage(:PerforatedCylinder,
  sysimage_path=joinpath(@__DIR__,"..","PerforatedCylinder_parallel.so"),
  precompile_execution_file=joinpath(@__DIR__,"warmup.jl"))
