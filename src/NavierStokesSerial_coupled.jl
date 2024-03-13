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
  println(dims)
  println(typeof(u1(0,0)))
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
  Y = MultiFieldFESpace([V, Q, Κ])
  Y₀ = MultiFieldFESpace([V, Q])

  # Define trial FESpaces from Dirichlet values
  U = TransientTrialFESpace(V, U0_dirichlet)
  P = TrialFESpace(Q)
  Η = TrialFESpace(Κ,[VectorValue(0.0,0.0),VectorValue(0.0,0.0),VectorValue(0.0,0.0)])
  X = TransientMultiFieldFESpace([U, P,Η])
  X₀ = MultiFieldFESpace([U(0.0), P])

  # Stokes for pre-initalize NS
  σ_dev_f(ε) = 2 * ν_f * ε #  Cauchy stress tensor for the fluid
  a((u, p), (v, q)) = ∫(ε(v) ⊙ (σ_dev_f ∘ ε(u)) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ_f
  l((v, q)) = ∫(0.0 * q)dΩ_f
  stokes_op = AffineFEOperator(a,l,X₀,Y₀)

  # Linear Solver
  to_logfile("Stokes solve")
  ls₀ = LUSolver()
  u_ST, p_ST = solve(ls₀,stokes_op)

  # initial condition NS
  to_logfile("Navier-Stokes operator")
  xh₀ = interpolate_everywhere([u_ST, p_ST, VectorValue(0.0,0.0)],X(0.0))
  vh₀ = interpolate_everywhere((u0(0),0.0,VectorValue(0.0,0.0)),X(0.0))

  # Explicit FE functions
  # global ηₙₕ = interpolate(u0(0),U(0.0))
  # global uₙₕ = interpolate(u_ST,U(0.0))
  # global fv_u = zero_free_values(U(0.0))

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  h = CellField(get_cell_measure(Ω_f),Ω_f)
  h2 = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω_f)),Ω_f)
  τₘ⁻¹(u) = (c₁*ν_f/h2 + c₂*((u⋅u).^(1/2))/h)
  τₘ(u) = 1/τₘ⁻¹(u)
  τc(u) = cc *(h2/(c₁*τₘ(u)))
  dτₘ(u,du) = -1.0/(τₘ⁻¹(u)*τₘ⁻¹(u)) * (c₂*(u⋅du)./((u⋅u).^(1/2)+1.0e-14))
  dτc(u,du) = -cc*h2/c₁ * (1/(τₘ(u)*τₘ(u))) * dτₘ(u,du)

  # Orthogonal projection
  # aη(u) = (η,κ) -> ∫( τₘ(u)*(η⋅κ) )dΩ_f
  # bη(u) = (κ) -> ∫( τₘ(u)*((∇(u)'⋅u)⋅κ) )dΩ_f
  # op_proj(u) = AffineFEOperator(aη(u),bη(u),Η,Κ)
  # ls_proj = LUSolver()
  # ηₕ(u) = solve(ls_proj,op_proj(u))

  # buffer = Ref{Any}((η=nothing,t=nothing))
  # function η(u,t)
  #   if buffer[].t == t
  #     return buffer[].η
  #   else
  #     buffer[] = (η=ηₕ(u),t=t)
  #     return buffer[].η
  #   end
  # end

  # Weak form
  c(a,u,v) = 0.5*((∇(u)'⋅a)⋅v - u⋅(∇(v)'⋅a))
  neg(a) = min(a,0.0)
  mass(t,(∂ₜu,),(v,)) = ∫( ∂ₜu⋅v )dΩ_f
  res(t,(u,p,η),(v,q,κ)) = ∫( c(u,u,v) )dΩ_f +
                          ∫( ε(v) ⊙ (σ_dev_f ∘ ε(u)) )dΩ_f -
                          ∫( p*(∇⋅v) )dΩ_f +
                          ∫( (∇⋅u)*q )dΩ_f +
                          ∫( τₘ(u)*((∇(u)'⋅u - η)⋅(∇(v)'⋅u-κ)) )dΩ_f +
                          ∫( τc(u)*((∇⋅u)*(∇⋅v)) )dΩ_f +
                       ∫( (u⋅v)*(0.5*(u⋅n_Γout)-neg∘(u⋅n_Γout)) )dΓout
  jac(t,(u,p,η),(du,dp,dη),(v,q,κ)) =
    ∫( c(du,u,v) )dΩ_f +
    ∫( c(u,du,v) )dΩ_f +
    ∫( ε(v) ⊙ (σ_dev_f ∘ ε(du)) )dΩ_f -
    ∫( dp*(∇⋅v) )dΩ_f +
    ∫( (∇⋅du)*q )dΩ_f +
    ∫( τₘ(u)*((∇(u)'⋅u - η)⋅(∇(v)'⋅du)) )dΩ_f +
    ∫( τₘ(u)*((∇(du)'⋅u + ∇(u)'⋅du - dη)⋅(∇(v)'⋅u-κ)) )dΩ_f +
    ∫( τc(u)*((∇⋅du)*(∇⋅v)) )dΩ_f +
    ∫( dτₘ(u,du)*((∇(u)'⋅u - η)⋅(∇(v)'⋅u-κ)) )dΩ_f +
    ∫( dτc(u,du)*((∇⋅u)*(∇⋅v)) )dΩ_f +
    ∫( (du⋅v)*(0.5*(u⋅n_Γout)-neg∘(u⋅n_Γout)) )dΓout +
    ∫( (u⋅v)*(0.5*(du⋅n_Γout)-neg∘(du⋅n_Γout)) )dΓout
  jac_t(t,(u,),(dut,),(v,)) = ∫( dut⋅v )dΩ_f

  # NS operator
  op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y;constant_mass=true)
  # op = TransientSemilinearFEOperator(mass, res, X, Y;constant_mass=true)

  # Nonlinear Solver
  nls = NLSolver(LUSolver(),show_trace=true,method=:newton,iterations=10,ftol=1.0e-6)
  ls_mass = LUSolver()

  # ODE solvers:
  # 1 time step with BE to kill spurious oscillations in force
  ode_solver₁ = ThetaMethod(nls,Δt,1.0)
  ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_Midpoint_1_2()))
  # ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_3_3()))
  # ode_solver₂ = EXRungeKutta(ls_mass,Δt,ButcherTableau(EXRK_RungeKutta_4_4()))

  xₜ₁ = solve(ode_solver₁,op,t₀,t₀+Δt,xh₀)
  function get_step(xₜ)
    for (t,x) in xₜ
      return x
    end
  end
  xh₁ = get_step(xₜ₁)
  xₜ = solve(ode_solver₂,op,t₀+Δt,tf,xh₁)

  # Postprocess
  println("Postprocess")
  global tout = 0
  createpvd("NS_test") do pvd
    for (t,(uh,ph,ηₕ)) in xₜ
      to_logfile("Time: $t")
      to_logfile("=======================")
      Fx, Fy = sum(∫((n_ΓS ⋅ σ_dev_f(ε(uh))) - ph * n_ΓS) * dΓₛ)
      to_forcefile(t,Fx,Fy)
      # uₙₕ = interpolate!(uh,fv_u,U(t))
      # ηₙₕ = solve(ls_proj,op_proj(uₙₕ))
      if t>=tout
        pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"eta_n"=>ηₕ],order=2)
        tout=t+Δtout
      end
    end
  end

  close(io)
  close(io_force)

  return nothing
end
