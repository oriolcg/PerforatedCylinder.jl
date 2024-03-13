function run_test_parallel(parts,mesh_file::String,force_file::String,Δt,tf,Δtout)

  if i_am_main(parts)
    io = open("output.log", "w")
    forces_path=ENV["PerforatedCylinder_FORCES"]
    full_force_path = joinpath(forces_path,force_file)
    io_force = open(full_force_path, "w")
  end
  function to_logfile(x...)
    if i_am_main(parts)
      write(io,join(x, " ")...)
      write(io,"\n")
      flush(io)
    end
  end
  function to_forcefile(x...)
    if i_am_main(parts)
      write(io_force,join(x, " ")...)
      write(io_force,"\n")
      flush(io_force)
    end
  end

  # Geometry
  to_logfile("Geometry")
  DIRICHLET_tags = ["inlet", "walls", "monopile"]
  DIRICHLET_masks = [(true,true),(false,true),(true,true)]
  FLUID_LABEL = "fluid"
  OUTLET_LABEL = "outlet"
  meshes_path=ENV["PerforatedCylinder_MESHES"]
  full_mesh_path = joinpath(meshes_path,mesh_file)
  to_logfile("Mesh file: ",full_mesh_path)
  testname = replace(mesh_file,".msh" =>"")
  model =  GmshDiscreteModel(parts,full_mesh_path)
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
  # D = 10 # m
  # Re = 1.0e6 #0
  H = 24 # m
  μ_f = 1.0#1.0e-3# rho * Vinf * D / Re #0.01 # Fluid viscosity
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
  Y = MultiFieldFESpace([V, Q])
  Κ = TestFESpace(Ω, reffeᵤ,  dirichlet_tags = DIRICHLET_tags, dirichlet_masks=DIRICHLET_masks, conformity = :H1)

  # Define trial FESpaces from Dirichlet values
  U = TransientTrialFESpace(V, U0_dirichlet)
  P = TrialFESpace(Q)
  X = TransientMultiFieldFESpace([U, P])
  Η = TrialFESpace(Κ,[u0(0.0), u0(0.0), u0(0.0)])

  # Stokes for pre-initalize NS
  σ_dev_f(ε) = 2 * ν_f * ε #  Cauchy stress tensor for the fluid
  a((u, p), (v, q)) = ∫(ε(v) ⊙ (σ_dev_f ∘ ε(u)) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ_f
  l((v, q)) = ∫(0.0 * q)dΩ_f
  stokes_op = AffineFEOperator(a,l,X(0.0),Y)

  # Setup solver via low level PETSC API calls
  function mykspsetup(ksp)
    pc       = Ref{GridapPETSc.PETSC.PC}()
    mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
    @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
    @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
    @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
    @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
    @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 2)
    @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  7, 0)
    @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  14, 500000)
    @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  24, 1)
    # @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
    # @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
    @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-10)
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
  end

  # Linear Solver
  to_logfile("Stokes solve")
  ls₀ = PETScLinearSolver(mykspsetup)
  u_ST, p_ST = solve(ls₀,stokes_op)
  # u_ST, p_ST = solve(stokes_op)

  # initial condition NS
  to_logfile("Navier-Stokes operator")
  xh₀ = interpolate_everywhere([u_ST, p_ST],X(0.0))
  vh₀ = interpolate_everywhere((u0(0),0.0),X(0.0))

  # Explicit FE functions
  global ηₙₕ = interpolate(u0(0),U(0.0))
  global uₙₕ = interpolate(u_ST,U(0.0))
  global fv_u = zero_free_values(U(0.0))

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  h2map = map(Ω_f.trians) do trian
    CellField(get_cell_measure(trian),trian)
  end
  h2 = DistributedCellField(h2map,Ω_f)
  hmap = map(Ω_f.trians) do trian
    CellField(lazy_map(dx->dx^(1/2),get_cell_measure(trian)),trian)
  end
  h = DistributedCellField(hmap,Ω_f)
  abs_(u) = (u⋅u).^(1/2) .+ 1e-14
  τₘ⁻¹(u) = (c₁*ν_f/h2 + c₂*(abs_(u))/h)
  τₘ(u) = 1.0 / τₘ⁻¹(u)
  τc(u) = cc *(h2/(c₁*τₘ(u)))

  dabs_(u,du) = (1/abs_(u))*(u⋅du)
  dτₘ⁻¹(u,du) = c₂/h*dabs_(u,du)
  dτₘ(u,du) = -1/(τₘ⁻¹(u)*τₘ⁻¹(u)) * dτₘ⁻¹(u,du)
  dτc(u,du) = -cc*h2/c₁ * (1/(τₘ(u)*τₘ(u))) * dτₘ(u,du)


  # Weak form
  c(a,u,v) = 0.5*((∇(u)'⋅a)⋅v - u⋅(∇(v)'⋅a))
  mass(t,(∂ₜu,),(v,)) = ∫( ∂ₜu⋅v )dΩ_f
  res(t,(u,p),(v,q)) = ∫( c(u,u,v) + ε(v) ⊙ (σ_dev_f ∘ ε(u)) - p*(∇⋅v) + (∇⋅u)*q +
                          τₘ(u)*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅u)) + τc(u)*((∇⋅u)*(∇⋅v)) )dΩ_f +
                       ∫( 0.5*(u⋅v)*(u⋅n_Γout) )dΓout
  jac(t,(u,p),(du,dp),(v,q)) = ∫( c(du,u,v) + c(u,du,v) + ε(v) ⊙ (σ_dev_f ∘ ε(du)) - dp*(∇⋅v) + (∇⋅du)*q +
                                  τₘ(u)*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅du) + (∇(du)'⋅u + ∇(u)'⋅du)⋅(∇(v)'⋅u)) +
                                  τc(u)*((∇⋅du)*(∇⋅v)) +
                                  dτₘ(u,du)*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅u)) +
                                  dτc(u,du)*((∇⋅u)*(∇⋅v)) )dΩ_f +
                                ∫( 0.5*((du⋅v)*(u⋅n_Γout)+(u⋅v)*(du⋅n_Γout)) )dΓout
  jac_t(t,(u,p),(dut,dpt),(v,q)) = ∫( dut⋅v )dΩ_f

  # Orthogonal projection
  aη(u) = (η,κ) -> ∫( τₘ(u)*(η⋅κ) )dΩ_f
  bη(u) = (κ) -> ∫( τₘ(u)*((∇(u)'⋅u)⋅κ) )dΩ_f
  op_proj(u) = AffineFEOperator(aη(u),bη(u),Η,Κ)
  ls_proj = PETScLinearSolver(mykspsetup)

  # NS operator
  op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y;constant_mass=true)
  # op = TransientSemilinearFEOperator(mass, res, X, Y;constant_mass=true)

  # Nonlinear Solver
  nls = PETScNonlinearSolver()
  ls_mass = PETScLinearSolver(mykspsetup)

  # Nonlinear Solver
  #nls = NLSolver(ls,show_trace=true,method=:newton,iterations=10)

  # ODE solvers:
  # 1 time step with BE to kill spurious oscillations in force
  ode_solver₁ = ThetaMethod(nls,Δt,1.0)
  ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_Midpoint_1_2()))
  # ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_3_3()))

  xₜ₁ = solve(ode_solver₁,op,t₀,t₀+Δt,xh₀)
  function get_step(xₜ)
    for (t,x) in xₜ
      return x
    end
  end
  xh₁ = get_step(xₜ₁)
  xₜ = solve(ode_solver₂,op,t₀+Δt,tf,xh₁)

  # Postprocess
  if i_am_main(parts)
    println("Postprocess")
  end
  global tout = 0
  createpvd(parts,"NS_test") do pvd
    for (t,(uh,ph)) in xₜ
      to_logfile("Time: $t")
      to_logfile("=======================")
      Fx, Fy = sum(∫((n_ΓS ⋅ σ_dev_f(ε(uh))) - ph * n_ΓS) * dΓₛ)
      to_forcefile(t,Fx,Fy)
      if t>tout
        pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"un"=>uₙₕ,"eta_n"=>ηₙₕ])
        tout=t+Δtout
      end
      uₙₕ = interpolate!(uh,fv_u,U(t))
      ηₙₕ = solve(ls_proj,op_proj(uₙₕ))
    end
  end

  if i_am_main(parts)
    close(io)
    close(io_force)
  end

  return nothing
end
