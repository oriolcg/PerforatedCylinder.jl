using PerforatedCylinder

println("generating meshes...")
nbeta = 21
nalpha = 1
nperfs = 40
perf_cases = 3:30
porosities = 0.3:0.02:0.7
alphas = [0.0]#:15.0/(nalpha-1):15.0
for num_perforations in perf_cases
  for β in porosities
    for α in alphas
      println(num_perforations,"-",β,"-",α)
      PerforatedCylinder.create_mesh(num_perforations = num_perforations, β = β, α = α)
    end
  end
end
