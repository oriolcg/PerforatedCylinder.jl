function create_mesh(;
  # Domain Parameters
  L = 10,
  H = 4,
  D = 1,
  t = 0.05,
  Cx = 4,
  Cy = 2,
  # Perforations
  num_perforations = 3,
  β = 0.3,
  α = 30,
  # Mesh Parameters
  h_coarse = 0.3,
  h_fine = 1.0e-2,
  dxLeft = 0.5,
  dxRight = 6,
  dyTop = 0.5,
  dyBottom = 0.5,
  decay_factor_left = 5,
  decay_factor_right = 1.5,
  decay_exponent_left = 0.8,
  decay_exponent_right = 0.8,
  )

  # Domain Parameters
  R = D/2+t/2
  r = D/2-t/2

  # Perforations
  ϕ = 2π/num_perforations
  θ = β*ϕ

  # Initialize
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)
  gmsh.model.add("t1")

  # monopile
  gmsh.model.occ.addDisk(Cx,Cy,0,R,R,2)
  gmsh.model.occ.addDisk(Cx,Cy,0,r,r,3)
  monopile = gmsh.model.occ.cut((2,2), (2,3), 4)

  # Perforations
  w = θ*R
  l = R+2t
  perforations = []
  for i in 1:num_perforations
    tmp1 = gmsh.model.occ.addRectangle(Cx-w/2,Cy-l,0,w,l)
    gmsh.model.occ.rotate((2,tmp1),Cx,Cy,0,0,0,1,(i-1)*ϕ+α/180*π)
    tmp2 = gmsh.model.occ.intersect((2,tmp1),monopile[1][end],-1,true,false)
    push!(perforations,tmp2)
  end
  global tmp3 = monopile
  for perforation in perforations
    tmp3 = gmsh.model.occ.cut(tmp3[1][end],perforation[1][end],-1,true,true)
  end
  monopile_pieces = gmsh.model.occ.getEntities(2)

  # Background Domain
  global domain = gmsh.model.occ.addRectangle(0,0,0,L,H)

  # Subtract pieces from domain
  for piece in monopile_pieces
    domain = gmsh.model.occ.cut((2,domain),piece,-1,true,true)
    domain = domain[1][end][end]
  end

  # Get entities
  fluid_domain = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,L+0.1,H+0.1,0.1,2)
  wall1 = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,H-0.1,-0.1,L+0.1,H+0.1,0.1,1)
  wall2 = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,L+0.1,0.1,0.1,1)
  inlet_point = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,0.1,H+0.1,0.1,0)
  inlet_line = gmsh.model.occ.getEntitiesInBoundingBox(-0.1,-0.1,-0.1,0.1,H+0.1,0.1,1)
  outlet_point = gmsh.model.occ.getEntitiesInBoundingBox(L-0.1,-0.1,-0.1,L+0.1,H+0.1,0.1,0)
  outlet_line = gmsh.model.occ.getEntitiesInBoundingBox(L-0.1,-0.1,-0.1,L+0.1,H+0.1,0.1,1)
  monopile_point = gmsh.model.occ.getEntitiesInBoundingBox(Cx-R-0.1,Cy-R-0.1,-0.1,Cx+R+0.1,Cy+R+0.1,0.1,0)
  monopile_line = gmsh.model.occ.getEntitiesInBoundingBox(Cx-R-0.1,Cy-R-0.1,-0.1,Cx+R+0.1,Cy+R+0.1,0.1,1)


  # Get entity tags
  fluid_domain_tags = [ entity[2] for entity in fluid_domain]
  wall_tags = [ entity[2] for entity in wall1]
  append!(wall_tags , [ entity[2] for entity in wall2])
  inlet_point_tags = [ entity[2] for entity in inlet_point]
  inlet_line_tags = [ entity[2] for entity in inlet_line]
  outlet_point_tags = [ entity[2] for entity in outlet_point]
  outlet_line_tags = [ entity[2] for entity in outlet_line]
  monopile_point_tags = [ entity[2] for entity in monopile_point]
  monopile_line_tags = [ entity[2] for entity in monopile_line]

  println("fluid_domain", fluid_domain)

  # Physical group
  gmsh.model.addPhysicalGroup(0,inlet_point_tags,1,"inlet")
  gmsh.model.addPhysicalGroup(1,inlet_line_tags,1,"inlet")
  gmsh.model.addPhysicalGroup(0,outlet_point_tags,2,"outlet")
  gmsh.model.addPhysicalGroup(1,outlet_line_tags,2,"outlet")
  gmsh.model.addPhysicalGroup(1,wall_tags,3,"walls")
  gmsh.model.addPhysicalGroup(0,monopile_point_tags,4,"monopile")
  gmsh.model.addPhysicalGroup(1,monopile_line_tags,4,"monopile")
  pg1 = gmsh.model.addPhysicalGroup(2,fluid_domain_tags)#,5,"fluid")
  gmsh.model.setPhysicalName(2,pg1,"fluid")

  # Synchronize
  gmsh.model.occ.synchronize()
  gmsh.model.geo.synchronize()

  # Define mesh size
  function meshSizeCallback(dim,tag,x,y,z,lc)
    if (Cx-R-dxLeft)<x<(Cx) && √((x-Cx)^2+(y-Cy)^2) < (R+dxLeft)
      dist = abs(√((x-Cx)^2+(y-Cy)^2) - R)/R
      return min(h_fine * (1 + decay_factor_left * (dist^decay_exponent_left)), h_coarse)
    elseif (Cx)<x<(Cx+R+dxRight) && (Cy-R-dyBottom)<y<(Cy+R+dyTop)
      dist = abs(√((x-Cx)^2+(y-Cy)^2) - R)/R
      return min(h_fine * (1 + decay_factor_right * (dist^decay_exponent_right)), h_coarse)
    else
      return h_coarse
    end
  end
  gmsh.model.mesh.setSizeCallback(meshSizeCallback)
  gmsh.model.mesh.generate()

  println(gmsh.model.getEntitiesForPhysicalGroup(2,5))

  # Finalize
  β2 = round(β;digits=2)
  α2 = round(α,digits=2)
  filename = "$num_perforations-$β2-$α2.msh"
  meshes_path=ENV["CNN_NS_MESHES"]
  mesh_file = joinpath(meshes_path,filename)
  gmsh.write(mesh_file)
  gmsh.finalize()

end
