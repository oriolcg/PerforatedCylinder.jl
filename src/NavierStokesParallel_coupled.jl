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
  Y = MultiFieldFESpace([V, Q, Κ])
  Y₀ = MultiFieldFESpace([V, Q])

  # Define trial FESpaces from Dirichlet values
  U = TransientTrialFESpace(V, U0_dirichlet)
  P = TrialFESpace(Q)
  Η = TrialFESpace(Κ,[VectorValue(0.0,0.0),VectorValue(0.0,0.0),VectorValue(0.0,0.0)])
  X = TransientMultiFieldFESpace([U, P,Η])
  X₀ = MultiFieldFESpace([U(0.0), P])

  # Stokes for pre-initalize NS
  a((u, p), (v, q)) = ∫( 2ν_f*(ε(v) ⊙ ε(u)) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ_f
  l((v, q)) = ∫(0.0 * q)dΩ_f
  stokes_op = AffineFEOperator(a,l,X₀,Y₀)

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

  # initial condition NS
  to_logfile("Navier-Stokes operator")
  xh₀ = interpolate_everywhere([u_ST, p_ST, VectorValue(0.0,0.0)],X(0.0))

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  τₘ, τc, dτₘ, dτc = get_stabilization_parameters(Ω_f, ν_f, c₁, c₂, cc)

  mass(t,(∂ₜu,),(v,)) = ∫( ∂ₜu⋅v )dΩ_f
  res(t,(u,p,η),(v,q,κ)) = conv(u,u,v,dΩ_f) +
                           lap(ν_f,u,v,dΩ_f) -
                           div(v,p,dΩ_f) +
                           div(u,q,dΩ_f) +
                           stab(τₘ,u,u,p,η,v,q,κ,dΩ_f) +
                           graddiv(τc,u,u,v,dΩ_f) +
                           cΓ(u,u,v,n_Γout,dΓout)
  jac(t,(u,p,η),(du,dp,dη),(v,q,κ)) =
    lap(ν_f,du,v,dΩ_f) -
    div(v,dp,dΩ_f) +
    div(du,q,dΩ_f) +
    dconv(u,u,du,du,v,dΩ_f) +
    dstab(τₘ,dτₘ,u,u,p,η,du,du,dp,dη,v,q,κ,dΩ_f) +
    dgraddiv(τc,dτc,u,u,du,du,v,dΩ_f) +
    ∫( (du⋅v)*(0.5*(u⋅n_Γout)-neg∘(u⋅n_Γout)) )dΓout +
    ∫( (u⋅v)*(0.5*(du⋅n_Γout)-neg∘(du⋅n_Γout)) )dΓout
  jac_t(t,(u,),(dut,),(v,)) = ∫( dut⋅v )dΩ_f

  # NS operator
  op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y;constant_mass=true)

  # Nonlinear Solver
  nls = PETScNonlinearSolver()
  ls_mass = PETScLinearSolver(mykspsetup)

  # ODE solvers:
  # 1 time step with BE to kill spurious oscillations in force
  ode_solver₁ = ThetaMethod(nls,Δt,1.0)
  ode_solver₂ = DIMRungeKutta(nls,ls_mass,Δt,ButcherTableau(SDIRK_Midpoint_1_2()))

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
    for (t,(uh,ph,ηₕ)) in xₜ
      to_logfile("Time: $t")
      to_logfile("=======================")
      Fx, Fy = sum(∫(2ν_f*(n_ΓS ⋅ ε(uh)) - ph * n_ΓS) * dΓₛ)
      to_forcefile(t,Fx,Fy)
      if t>=tout
        pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"eta_n"=>ηₕ,"usgs"=>uₛ(τₘ,uh,uh,ph,ηₕ)],order=2)
        tout=t+Δtout
      end
    end
  end

  if i_am_main(parts)
    close(io)
    close(io_force)
  end

  return nothing
end
