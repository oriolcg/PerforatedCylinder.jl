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
  DIRICHLET_tags = ["inlet", "walls"]#, "monopile"]
  DIRICHLET_masks = [(true,true),(false,true)]#,(true,true)]
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
  U0_dirichlet = [u1, u1]#, u0]
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
  Η = TrialFESpace(Κ,[VectorValue(0.0,0.0),VectorValue(0.0,0.0)])#,VectorValue(0.0,0.0)])
  X = TransientMultiFieldFESpace([U, P])
  X₀ = MultiFieldFESpace([U(0.0), P])

  # Stokes for pre-initalize NS
  a((u, p), (v, q)) = ∫( 2ν_f*(ε(v)⊙ε(u)) )dΩ_f - ∫( (∇⋅v)*p + q*(∇⋅u))dΩ_f
  l((v, q)) = 0
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
  h, h2 = get_mesh_sizes(Ω)
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
      println("Computing uₙ")
      buffer2[] = (u=uₙₕ(u,t),t=t)
      return buffer2[].u
    end
  end

  # Wall law
  χ = 0.4
  B = 5.5
  γ = 100.0
  Cᵦ = 32.0
  fₜ(u) = u + exp(-χ*B)*(exp(χ*u) - 1.0 - χ*u - 0.5*(χ^2*u^2) -(1/6)*(χ^3*u^3) )
  dfₜ(u) = 1 + χ*exp(-χ*B)*(exp(χ*u) - 1.0 - χ*u - 0.5*(χ^2*u^2) )
  τ_(h) = Cᵦ*ν_f/h
  y⁺_(h,absᵤ,τ) = h/(γ*ν_f)*√(absᵤ*τ)
  I = TensorValue(1.0,0.0,0.0,1.0)
  function τᵦ(u,u₀,h,n)
    uₜ = (I-n⊗n)⋅(u-u₀)
    absᵤ = norm(uₜ)
    τ = τ_(h)
    u⁺ = √(absᵤ/τ)
    y⁺ = u⁺
    r = y⁺ - fₜ(u⁺)
    while abs(r) > 1.0e-6
      dr = h/(2ν_f*γ)*√(absᵤ/τ) + dfₜ(u⁺)*0.5*(√(absᵤ/(τ^3)))
      dτ = -r/dr
      τ += dτ
      u⁺ = √(absᵤ/τ)
      y⁺ = y⁺_(h,absᵤ,τ)
      r = y⁺ - fₜ(u⁺)
    end
    τₜ = τ*(I-n⊗n)
    τₙ = γ*τ_(h)*(n⊗n)
    return τₜ+τₙ
  end
  buffer3 = Ref{Any}((τᵦ=nothing,t=nothing))
  function τᵦ(u,u₀,h,n,t)
    if buffer3[].t == t
      return buffer3[].τᵦ
    else
      println("Computing τᵦ")
      buffer3[] = (τᵦ=τᵦ∘(u,u₀,h,n),t=t)
      return buffer3[].τᵦ
    end
  end

  # TMP weak form
  # stab_expl_(u,p,v,q,t) = ∫( τₘ(u,h,h2) * ((∇(u)'⋅u + ∇(p) - η(uₙ(u,t),t))⋅(∇(v)'⋅u)+∇(q)) )dΩ_f
  𝒞ᵤ(a,∇u) = a⋅∇u
  ℒᵤ(a,∇u,∇p) = 𝒞ᵤ(a,∇u) + ∇p
  𝒫ᵤ(a,∇u,∇p,η) = ℒᵤ(a,∇u,∇p) - η
  skew_conv(a,u,v,∇u,∇v) =  0.5*(𝒞ᵤ(a,∇u)⋅v - 𝒞ᵤ(a,∇v)⋅u)
  sym_lapl(εu,εv) = 2ν_f*(εu⊙εv)
  div_term(divu,q) = divu*q
  skew_conv_Γ(a,u,v,n) = (u⋅v)*(0.5*(a⋅n)-neg(a⋅n))
  penalty(τ,u₀,u,v) = (τ⋅(u-u₀))⋅v
  dpenalty(τ,du,v) = (τ⋅du)⋅v
  complementary_uv(u₀,u,εv,n) = 2ν_f * ((n⋅εv)⋅(u-u₀))
  complementary_uq(u₀,u,q,n) = 2ν_f * ((q*n)⋅(u-u₀))
  dcomplementary_uv(du,εv,n) = 2ν_f * ((n⋅εv)⋅(du))
  dcomplementary_uq(du,q,n) = 2ν_f * ((q*n)⋅(du))
  u0cf(t) = CellField(u0(t),Γ_S)

  # Weak form
  mass(t,(∂ₜu,),(v,)) = ∫( ∂ₜu⋅v )dΩ_f
  res(t,(u,p),(v,q)) =
    ∫( (sym_lapl∘(ε(u),ε(v))) +
       (skew_conv∘(u,u,v,∇(u),∇(v))) +
       (div_term∘((∇⋅u),q)) -
       (div_term∘((∇⋅v),p)) +
       (τₘ∘(uₙ(u,t),h,h2)) * ((𝒫ᵤ∘(uₙ(u,t),∇(u),∇(p),η(uₙ(u,t),t)))⋅((𝒞ᵤ∘(uₙ(u,t),∇(v)))+∇(q))) +
       (τc∘(uₙ(u,t),h,h2)) * ((∇⋅u)*(∇⋅v)) )dΩ_f +
    ∫( (skew_conv_Γ∘(u,u,v,n_Γout)) )dΓout +
    ∫( (penalty∘(τᵦ(uₙ(u,t),u0cf(t),h,n_ΓS,t),u0cf(t),u,v)) -
       (complementary_uv∘(u0(t),u,ε(v),n_ΓS)) +
       (complementary_uq∘(u0(t),u,q,n_ΓS)) -
       (dcomplementary_uv∘(v,ε(u),n_ΓS)) +
       (dcomplementary_uq∘(v,p,n_ΓS)) )dΓₛ
  jac(t,(u,p),(du,dp),(v,q)) =
    ∫( (sym_lapl∘(ε(du),ε(v)))  +
        0.5*((𝒞ᵤ∘(du,∇(u)))⋅v + (𝒞ᵤ∘(u,∇(du)))⋅v - (𝒞ᵤ∘(du,∇(v)))⋅u - (𝒞ᵤ∘(u,∇(v)))⋅du) +
       (div_term∘((∇⋅du),q)) -
       (div_term∘((∇⋅v),dp)) +
       (τₘ∘(uₙ(u,t),h,h2)) * ((𝒞ᵤ∘(uₙ(u,t),∇(du)))+∇(dp))⋅((𝒞ᵤ∘(uₙ(u,t),∇(v)))+∇(q)) +
       (τc∘(uₙ(u,t),h,h2)) * ((∇⋅du)*(∇⋅v)) )dΩ_f +
    ∫( (du⋅v)*(0.5*(u⋅n_Γout)-neg∘(u⋅n_Γout)) )dΓout +
    ∫( (u⋅v)*(0.5*(du⋅n_Γout)-neg∘(du⋅n_Γout)) )dΓout +
    ∫( (dpenalty(τᵦ(uₙ(u,t),u0cf(t),h,n_ΓS,t),du,v)) )dΓₛ -
    ∫( (dcomplementary_uv(du,ε(v),n_ΓS)) )dΓₛ +
    ∫( (dcomplementary_uq(du,q,n_ΓS)) )dΓₛ -
    ∫( (dcomplementary_uv(v,ε(du),n_ΓS)) )dΓₛ +
    ∫( (dcomplementary_uq(v,dp,n_ΓS)) )dΓₛ
  jac_t(t,(u,),(dut,),(v,)) = ∫( dut⋅v )dΩ_f

  # NS operator
  op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y;constant_mass=true)

  # Nonlinear Solver
  nls = NLSolver(LUSolver(),show_trace=true,method=:newton,iterations=10,ftol=1.0e-6)#, linesearch=BackTracking())
  ls = LUSolver()

  # ########################################
  # # IMEX operators
  # stiffness(t,(u,p),(v,q)) =
  #   ∫( (sym_lapl∘(ε(u),ε(v))) +
  #      (div_term∘((∇⋅u),q)) -
  #      (div_term∘((∇⋅v),p)) )dΩ_f
  # im_jac(t,(u,p),(du,dp),(v,q)) = stiffness(t,(du,dp),(v,q))
  # ex_res(t,(u,p),(v,q)) =
  #   ∫( (skew_conv∘(u,u,v,∇(u),∇(v))) +
  #       (τₘ∘(u,h,h2)) * ((𝒫ᵤ∘(u,∇(u),∇(p),η(u,t)))⋅((𝒞ᵤ∘(u,∇(v)))+∇(q))) +
  #       (τc∘(u,h,h2)) * ((∇⋅u)*(∇⋅v)) )dΩ_f +
  #   ∫( (skew_conv_Γ∘(u,u,v,n_Γout)) )dΓout
  # ex_jac(t,(u,p),(du,dp),(v,q)) = ∫(0.0*(du⋅v))dΩ_f
  # im_op = TransientLinearFEOperator((stiffness, mass), (t,y)->0, (im_jac, jac_t), X, Y)#; constant_forms=(true,true))
  # ex_op = TransientFEOperator(ex_res, (ex_jac,), X, Y;assembler=Gridap.ODEs.get_assembler(im_op))
  # println(Gridap.ODEs.get_assembler(im_op))
  # println(Gridap.ODEs.get_assembler(ex_op))
  # imex_op = TransientIMEXFEOperator(im_op, ex_op)
  # ########################################

  # ODE solvers:
  # 1 time step with BE to kill spurious oscillations in force
  ode_solver₁ = ThetaMethod(nls,Δt,1.0)
  # ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_Midpoint_1_2()))
  ode_solver₂ = GeneralizedAlpha1(nls,Δt,0.5)
  # ode_solver₃ = RungeKutta(ls,ls,Δt,:IMEXRK_2_3_2)

  xₜ₁ = solve(ode_solver₁,op,t₀,t₀+Δt,xh₀)
  function get_step(xₜ)
    for (t,x) in xₜ
      return x
    end
  end
  xh₁ = get_step(xₜ₁)
  xh₁ₜ = (xh₁,interpolate_everywhere([VectorValue(0.0,0.0),0.0],X(t₀+Δt)))
  xₜ = solve(ode_solver₂,op,t₀+Δt,tf,xh₁ₜ)
  # xₜ = solve(ode_solver₃,imex_op,t₀+Δt,tf,xh₁)

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
