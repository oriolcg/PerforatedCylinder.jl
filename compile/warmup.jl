using PerforatedCylinder
PerforatedCylinder.main_serial()
# PerforatedCylinder.generate_meshes()
PerforatedCylinder.main_parallel(1)

# testname = "3D_monopile_coarse"
# mesh_file = "3D_monopile_coarse.msh"
# force_file = "3D_monopile_coarse.csv"
# PerforatedCylinder.main_parallel(1;
#   mesh_file=mesh_file,
#   force_file=force_file,
#   Δt=0.05,
#   tf=0.05,
#   Δtout=0.05,
# )
