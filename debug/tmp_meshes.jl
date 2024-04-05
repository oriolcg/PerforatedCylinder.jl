module TMP_meshes
using PerforatedCylinder

ENV["PerforatedCylinder_MESHES"] = "./debug/tmp_meshes/"
# mesh1 = PerforatedCylinder.create_mesh(
#   # Domain Parameters
#   L = 18,
#   H = 7,
#   D = 1,
#   t = 0.05,
#   Cx = 4,
#   Cy = 7/2,
#   # Perforations
#   num_perforations = 3,
#   β = 0.3,
#   α = 0,
#   # Mesh Parameters
#   h_coarse = 0.25,
#   h_fine = 7.5e-3,
#   dxLeft = 2,
#   dxRight = 18,
#   dyTop = 1.5,
#   dyBottom = 1.5,
#   decay_factor_left = 8,
#   decay_factor_right = 2,
#   decay_exponent_left = 0.8,
#   decay_exponent_right = 0.8
# )

# mesh2 = PerforatedCylinder.create_mesh(
#   # Domain Parameters
#   L = 18,
#   H = 7,
#   D = 1,
#   t = 0.05,
#   Cx = 4,
#   Cy = 7/2,
#   # Perforations
#   num_perforations = 3,
#   β = 0.3,
#   α = 0,
#   # Mesh Parameters
#   h_coarse = 0.3,
#   h_fine = 10.0e-3,
#   dxLeft = 2,
#   dxRight = 18,
#   dyTop = 1.5,
#   dyBottom = 1.5,
#   decay_factor_left = 8,
#   decay_factor_right = 2,
#   decay_exponent_left = 0.8,
#   decay_exponent_right = 0.8
# )

# mesh3 = PerforatedCylinder.create_mesh(
#   # Domain Parameters
#   L = 18,
#   H = 7,
#   D = 1,
#   t = 0.05,
#   Cx = 4,
#   Cy = 7/2,
#   # Perforations
#   num_perforations = 3,
#   β = 0.3,
#   α = 0,
#   # Mesh Parameters
#   h_coarse = 0.4,
#   h_fine = 15.0e-3,
#   dxLeft = 2,
#   dxRight = 18,
#   dyTop = 1.5,
#   dyBottom = 1.5,
#   decay_factor_left = 8,
#   decay_factor_right = 2,
#   decay_exponent_left = 0.8,
#   decay_exponent_right = 0.8
# )

tmp_coarse = PerforatedCylinder.create_mesh(
  # Domain Parameters
  L = 18,
  H = 7,
  D = 1,
  t = 0.05,
  Cx = 4,
  Cy = 7/2,
  # Perforations
  num_perforations = 3,
  β = 0.3,
  α = 0,
  # Mesh Parameters
  h_coarse = 2.0,
  h_fine = 1.0e-1,
  dxLeft = 2,
  dxRight = 18,
  dyTop = 1.5,
  dyBottom = 1.5,
  decay_factor_left = 8,
  decay_factor_right = 2,
  decay_exponent_left = 0.8,
  decay_exponent_right = 0.8
)


end
