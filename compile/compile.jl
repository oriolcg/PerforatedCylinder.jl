using PackageCompiler
create_sysimage(:CNN_NS,
  sysimage_path=joinpath(@__DIR__,"..","CNN_NS.so"),
  precompile_execution_file=joinpath(@__DIR__,"warmup.jl"))
