function run_test_parallel(parts,mesh_file::String,force_file::String,output_path,Δt,tf,Δtout)

  if i_am_main(parts)
    isdir(output_path) || mkdir(output_path)
    io = open(output_path*"/output.log", "w")
    io_force = open(force_file, "w")
  end
  function to_logfile(x...)
    if i_am_main(parts)
      write(io,join(x, " ")...)
      write(io,"\n")
      flush(io)
    end
  end

  # Geometry
  to_logfile("Geometry")
  DIRICHLET_tags = ["inlet", "walls", "monopile"]
  FLUID_LABEL = "fluid"
  OUTLET_LABEL = "outlet"
  meshes_path=ENV["CNN_NS_MESHES"]
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
  rho = 1.025 * 10^3 # kg/m^3
  Vinf = 1 # m/s
  D = 10 # m
  Re = 1.0e6 #0
  H = 24 # m
  μ_f = 1.0e-3# rho * Vinf * D / Re #0.01 # Fluid viscosity
  ν_f = μ_f / rho # kinematic viscosity

  # Boundary conditions and external loads
  u0(x, t) = VectorValue(0.0, 0.0)
  u0(t::Real) = x -> u0(x,t)
  #u1(x,t) = VectorValue( 1.5 * Vinf * x[2] * ( H - x[2] ) / ( (H/2)^2 ), 0.0 )
  u1(x,t) = VectorValue( Vinf, 0.0 )
  u1(t::Real) = x -> u1(x,t)
  U0_dirichlet = [u1, u1, u0]
  f(x) = VectorValue(0.0, 0.0)
  g(x) = 0.0

  # ODE solver
  t₀ = 0.0 # start [s]
  ρ∞ = 0.5

  to_logfile("FE spaces")
  # ReferenceFE
  reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)#,space=:P)
  reffeₚ = ReferenceFE(lagrangian, Float64, order - 1)#,space=:P)

  # Define test FESpaces
  V = TestFESpace(Ω_f, reffeᵤ,  dirichlet_tags = DIRICHLET_tags, conformity = :H1)
  Q = TestFESpace(Ω_f, reffeₚ,   conformity= :C0)
  Y = MultiFieldFESpace([V, Q])

  # Define trial FESpaces from Dirichlet values
  U = TransientTrialFESpace(V, U0_dirichlet)
  P = TrialFESpace(Q)
  X = TransientMultiFieldFESpace([U, P])

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
    @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  14, 5000)
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
  vh₀ = interpolate_everywhere((VectorValue(0.0,0.0),0.0),X(0.0))

  # Explicit FE functions
  global ηₙₕ = interpolate(VectorValue(0.0,0.0),U(0.0))
  global uₙₕ = interpolate(u_ST,U(0.0))
  global fv_u = zero_free_values(U(0.0))

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  # hmap = map_parts(Ω_f.trians) do trian
  #   get_cell_measure(trian)
  # end
  # h2 = DistributedCellField(hmap)
  # h = lazy_map(dx->dx^(1/2),h2)
  h = 0.05
  τₘ = 1/(c₁*ν_f/h^2 + c₂*(meas∘uₙₕ)/h)
  τc = cc *(h^2/(c₁*τₘ))

  # Weak form
  c(a,u,v) = 0.5*((∇(u)'⋅a)⋅v - u⋅(∇(v)'⋅a))
  res(t,(u,p),(v,q)) = ∫( ∂t(u)⋅v  + c(u,u,v) + ν_f*(∇(u)⊙∇(v)) - p*(∇⋅v) + (∇⋅u)*q +
                          τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅u)) + τc*((∇⋅u)*(∇⋅v)) )dΩ_f +
                       ∫( 0.5*(u⋅v)*(u⋅n_Γout) )dΓout
  jac(t,(u,p),(du,dp),(v,q)) = ∫( c(du,u,v) + c(u,du,v) + ν_f*(∇(du)⊙∇(v)) - dp*(∇⋅v) + (∇⋅du)*q +
                                  τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅du) + (∇(du)'⋅u + ∇(u)'⋅du)⋅(∇(v)'⋅u)) +
                                  τc*((∇⋅du)*(∇⋅v)) )dΩ_f +
                               ∫( 0.5*((du⋅v)*(u⋅n_Γout)+(u⋅v)*(du⋅n_Γout)) )dΓout
  jac_t(t,(u,p),(dut,dpt),(v,q)) = ∫( dut⋅v )dΩ_f
  op = TransientFEOperator(res,jac,jac_t,X,Y)

  # Orthogonal projection
  aη(η,κ) = ∫( τₘ*(η⋅κ) )dΩ_f
  bη(κ) = ∫( τₘ*((∇(uₙₕ)'⋅uₙₕ)⋅κ) )dΩ_f
  op_proj = AffineFEOperator(aη,bη,U(0.0),V)
  ls_proj = PETScLinearSolver()

  # # Define residual
  # res(t,(u, p), (v, q)) = ∫( ∂t(u)⋅v   + (∇(u)'⋅u) ⋅ v + ε(v) ⊙ (σ_dev_f ∘ ε(u))
  # - (∇ ⋅ v) * p + q * (∇ ⋅ u)   )dΩ_f#+ stab(u,p,v,q) )dΩ_f
  # jac(t,(u,p),(du,dp),(v,q)) = ∫( (∇(du)'⋅u) ⋅ v + (∇(u)'⋅du) ⋅ v + ε(v) ⊙ (σ_dev_f ∘ ε(du))
  # - (∇ ⋅ v) * dp + q * (∇ ⋅ du)   )dΩ_f
  # jac_t(t,(u, p),(dut,), (v, q)) = ∫( dut⋅v )dΩ_f

  # NS operator
  op = TransientFEOperator(res,jac,jac_t, X, Y)

  # Nonlinear Solver
  nls = PETScNonlinearSolver()

  # Nonlinear Solver
  #nls = NLSolver(ls,show_trace=true,method=:newton,iterations=10)

  # ODE solver
  ode_solver = GeneralizedAlpha(nls,Δt,ρ∞)

  xₜ = solve(ode_solver,op,(xh₀,vh₀),t₀,tf)

  # Postprocess
  if i_am_main(parts)
    println("Postprocess")
  end
  FD = Float64[]
  FL = Float64[]
  global tout = 0
   createpvd(parts,"NS_test") do pvd
    for ((uh,ph),t) in xₜ
      to_logfile("Time: $t")
      to_logfile("=======================")
      FR = sum(∫((n_ΓS ⋅ σ_dev_f(ε(uh))) - ph * n_ΓS) * dΓₛ)
      push!(FD,FR[1])
      push!(FL,FR[2])
      if t>tout
        pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"un"=>uₙₕ,"eta_n"=>ηₙₕ])
        tout=t+Δtout
      end
      uₙₕ = interpolate!(uh,fv_u,U(t))
      ηₙₕ = solve(ls_proj,op_proj)
    end
  end

  return FD, FL
end
