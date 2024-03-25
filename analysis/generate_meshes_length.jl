using PerforatedCylinder

println("generating meshes...")
perf_cases = [3,12,27]
porosities = [0.3,0.5,0.7]
lengths = 5:15
alphas = [0.0]#:15.0/(nalpha-1):15.0
for L in lengths
  for num_perforations in perf_cases
    for β in porosities
      α = 0.0
      println(L,"-",num_perforations,"-",β,"-",α)
      PerforatedCylinder.create_mesh_length(num_perforations = num_perforations, β = β, α = α, L=L)
    end
  end
end
