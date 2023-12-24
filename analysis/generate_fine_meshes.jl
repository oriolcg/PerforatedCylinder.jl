using PerforatedCylinder

println("generating meshes...")
nbeta = 21
nalpha = 1
nperfs = 40
perf_cases = 3:30
porosities = 0.3:0.02:0.7
alphas = [0.0]#:15.0/(nalpha-1):15.0

ENV["PerforatedCylinder_MESHES"] = "./data/fine_meshes/"
for num_perforations in perf_cases
  for β in porosities
    for α in alphas
      println(num_perforations,"-",β,"-",α)
      PerforatedCylinder.create_mesh(
        L = 7, Cx = 2, t = 0.1,
        num_perforations = num_perforations, β = β, α = α,
        h_coarse=0.4, h_fine=0.01, dxRight=5, dxLeft=0.5,
        decay_factor_right=10, decay_exponent_right=0.3,
        decay_factor_left=10, decay_exponent_left=0.6)
    end
  end
end
