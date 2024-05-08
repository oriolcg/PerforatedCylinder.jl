function run_test_serial(mesh_file::String,force_file::String,Δt,tf,Δtout)

  io = open("output.log", "w")
  forces_path=ENV["PerforatedCylinder_FORCES"]
  full_force_path = joinpath(forces_path,force_file)
  io_force = open(full_force_path, "w")
  function to_logfile(x...)
    write(io,join(x, " ")...)
    write(io,"\n")
    flush(io)
  end
  function to_forcefile(x...)
    write(io_force,join(x, " ")...)
    write(io_force,"\n")
    flush(io_force)
  end

  # Geometry
  to_logfile("Geometry")
  DIRICHLET_tags = ["inlet", "walls", "monopile"]
  DIRICHLET_masks = [(true,true),(false,true),(true,true)]
  meshes_path=ENV["PerforatedCylinder_MESHES"]
  full_mesh_path = joinpath(meshes_path,mesh_file)
  to_logfile("Mesh file: ",full_mesh_path)
  testname = replace(mesh_file,".msh" =>"")
  model =  GmshDiscreteModel(full_mesh_path)
  Ω = Triangulation(model)
  Ω_f = Triangulation(model, tags = "fluid")
  Γ_S = Boundary(model, tags = "monopile")
  Γ_out = Boundary(model, tags = "outlet")
  writevtk(model,testname)

  to_logfile("Measures")
  order = 2
  degree = 2 * order
  dΩ_f = Measure(Ω_f, degree)
  dΓₛ = Measure(Γ_S, degree)
  dΓout = Measure(Γ_out, degree)
  n_ΓS = get_normal_vector(Γ_S)
  n_Γout = get_normal_vector(Γ_out)

  # Physics parameters
  to_logfile("Parameters")
  rho = 1.0e3#1.025e3 # kg/m^3
  Vinf = 1 # m/s
  μ_f = 1.0e0# rho * Vinf * D / Re #0.01 # Fluid viscosity
  ν_f = μ_f / rho # kinematic viscosity

  # Boundary conditions and external loads
  dims = num_cell_dims(model)
  u0(x,t,::Val{2}) = VectorValue(0.0, 0.0)
  u1(x,t,::Val{2}) = VectorValue( Vinf, 0.0 )
  u0(x,t,::Val{3}) = VectorValue(0.0, 0.0, 0.0)
  u1(x,t,::Val{3}) = VectorValue( Vinf, 0.0, 0.0 )
  u0(x,t::Real) = u0(x,t,Val(dims))
  u1(x,t::Real) = u1(x,t,Val(dims))
  u0(t::Real) = x -> u0(x,t,Val(dims))
  u1(t::Real) = x -> u1(x,t,Val(dims))
  U0_dirichlet = [u1, u1, u0]
  g(x) = 0.0

  # ODE solver
  t₀ = 0.0 # start [s]
  ρ∞ = 0.5

  to_logfile("FE spaces")
  # ReferenceFE
  reffeᵤ = ReferenceFE(lagrangian, VectorValue{dims,Float64}, order)#,space=:P)
  reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)#,space=:P)

  # Define test FESpaces
  V = TestFESpace(Ω, reffeᵤ,  dirichlet_tags = DIRICHLET_tags, dirichlet_masks=DIRICHLET_masks, conformity = :H1)
  Q = TestFESpace(Ω, reffeₚ,   conformity= :C0)
  Κ = TestFESpace(Ω, reffeᵤ,  dirichlet_tags = DIRICHLET_tags, dirichlet_masks=DIRICHLET_masks, conformity = :H1)
  Y = MultiFieldFESpace([V, Q])
  # Y₀ = MultiFieldFESpace([V, Q])

  # Define trial FESpaces from Dirichlet values
  U = TransientTrialFESpace(V, U0_dirichlet)
  P = TrialFESpace(Q)
  Η = TrialFESpace(Κ,[VectorValue(0.0,0.0),VectorValue(0.0,0.0),VectorValue(0.0,0.0)])
  X = TransientMultiFieldFESpace([U, P])
  X₀ = MultiFieldFESpace([U(0.0), P])

  # Stokes for pre-initalize NS
  a((u, p), (v, q)) = ∫( 2ν_f*(ε(v) ⊙ ε(u)) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ_f
  l((v, q)) = ∫(0.0 * q)dΩ_f
  stokes_op = AffineFEOperator(a,l,X₀,Y)

  # Linear Solver
  to_logfile("Stokes solve")
  ls₀ = LUSolver()
  u_ST, p_ST = solve(ls₀,stokes_op)

  # initial condition NS
  to_logfile("Navier-Stokes operator")
  xh₀ = interpolate_everywhere([u_ST, p_ST],X(0.0))

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  h, h2 = get_mesh_sizes(Ω_f)
  τₘ, τc, dτₘ, dτc = get_stabilization_parameters_(ν_f, c₁, c₂, cc)

  # Orthogonal projection
  aη(u) = (η,κ) -> ∫( τₘ(u,h,h2)*(η⋅κ) )dΩ_f
  bη(u) = (κ) -> ∫( τₘ(u,h,h2)*((u⋅∇(u))⋅κ) )dΩ_f
  op_proj(u) = AffineFEOperator(aη(u),bη(u),Η,Κ)
  ls_proj = LUSolver()
  ηₕ(u) = solve(ls_proj,op_proj(u))
  fv_u = zero_free_values(U(0.0))
  uₙₕ(u,t) = interpolate!(u,fv_u,U(t-Δt))

  buffer = Ref{Any}((η=nothing,t=nothing))
  function η(u,t)
    if buffer[].t == t
      return buffer[].η
    else
      println("Computing η")
      buffer[] = (η=ηₕ(u),t=t)
      return buffer[].η
    end
  end
  buffer2 = Ref{Any}((u=nothing,t=nothing))
  function uₙ(u,t)
    if buffer2[].t == t
      return buffer2[].u
    else
      buffer2[] = (u=uₙₕ(u,t),t=t)
      println("Computing uₙ")
      return buffer2[].u
    end
  end

  # Weak form
  mass(t,(∂ₜu,),(v,)) = ∫( ∂ₜu⋅v )dΩ_f
  res(t,(u,p),(v,q)) =  conv(u,u,v,dΩ_f) +
                        lap(ν_f,u,v,dΩ_f) -
                        div(v,p,dΩ_f) +
                        div(u,q,dΩ_f) +
                        stab_expl(τₘ,h,h2,uₙ(u,t),u,p,η(uₙ(u,t),t),v,q,dΩ_f) +
                        graddiv(τc,h,h2,uₙ(u,t),u,v,dΩ_f) +
                        cΓ(u,u,v,n_Γout,dΓout)
  jac(t,(u,p),(du,dp),(v,q)) =
    lap(ν_f,du,v,dΩ_f) -
    div(v,dp,dΩ_f) +
    div(du,q,dΩ_f) +
    dconv(u,u,du,du,v,dΩ_f) +
    dstab_expl(τₘ,h,h2,uₙ(u,t),u,p,η(uₙ(u,t),t),du,dp,v,q,dΩ_f) +
    dgraddiv_expl(τc,h,h2,uₙ(u,t),u,du,v,dΩ_f) +
    ∫( (du⋅v)*(0.5*(u⋅n_Γout)-neg∘(u⋅n_Γout)) )dΓout +
    ∫( (u⋅v)*(0.5*(du⋅n_Γout)-neg∘(du⋅n_Γout)) )dΓout
  jac_t(t,(u,),(dut,),(v,)) = ∫( dut⋅v )dΩ_f

  # NS operator
  op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y;constant_mass=true)

  # Nonlinear Solver
  nls = NLSolver(LUSolver(),show_trace=true,method=:newton,iterations=10,ftol=1.0e-6)#, linesearch=BackTracking())
  ls_mass = LUSolver()

  # ODE solvers:
  # 1 time step with BE to kill spurious oscillations in force
  ode_solver₁ = ThetaMethod(nls,Δt,1.0)
  # ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_Midpoint_1_2()))
  ode_solver₂ = GeneralizedAlpha1(nls,Δt,0.5)

  xₜ₁ = solve(ode_solver₁,op,t₀,t₀+Δt,xh₀)
  function get_step(xₜ)
    for (t,x) in xₜ
      return x
    end
  end
  xh₁ = get_step(xₜ₁)
  xh₁ₜ = (xh₁,interpolate_everywhere([VectorValue(0.0,0.0),0.0],X(t₀+Δt)))
  xₜ = solve(ode_solver₂,op,t₀+Δt,tf,xh₁ₜ)

  # Postprocess
  println("Postprocess")
  global tout = 0
  createpvd("NS_test") do pvd
    for (t,(uh,ph)) in xₜ
      to_logfile("Time: $t")
      to_logfile("=======================")
      Fx, Fy = sum(∫(2ν_f*(n_ΓS ⋅ ε(uh)) - ph * n_ΓS) * dΓₛ)
      to_forcefile(t,Fx,Fy)
      ηh = buffer[].η
      if t>=tout
        pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"eta_n"=>ηh,"usgs"=>uₛ(τₘ,h,h2,uh,uh,ph,ηh)],order=2)
        tout=t+Δtout
      end
    end
  end

  close(io)
  close(io_force)

  return nothing
end
