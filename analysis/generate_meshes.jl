using CNN_NS

println("generating meshes...")
nbeta = 3
nalpha = 3
perf_cases = [3,9,27]
porosities = 0.3:(0.7-0.3)/(nbeta-1):0.7
alphas = 0.0:15.0/(nalpha-1):15.0
for num_perforations in perf_cases
  for β in porosities
    for α in alphas
      println(num_perforations,"-",β,"-",α)
      CNN_NS.create_mesh(num_perforations = num_perforations, β = β, α = α)
    end
  end
end
