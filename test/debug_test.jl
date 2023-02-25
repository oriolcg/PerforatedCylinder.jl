module debug_test

using Gridap
using Gridap.FESpaces: zero_free_values, interpolate!
using GridapGmsh
using Gridap.Fields: meas

full_mesh_path = joinpath(@__DIR__,"../data/meshes/tmp_coarse.msh")
full_force_path = joinpath(@__DIR__,"../data/forces/force.csv")
io = open("output.log", "w")
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
FLUID_LABEL = "fluid"
OUTLET_LABEL = "outlet"
to_logfile("Mesh file: ",full_mesh_path)
testname = "tmp_coarse"
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
rho = 1.025e3 # kg/m^3
Vinf = 1 # m/s
D = 10 # m
# Re = 1.0e6 #0
H = 24 # m
μ_f = 1.0e-3# rho * Vinf * D / Re #0.01 # Fluid viscosity
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
V = TestFESpace(Ω, reffeᵤ,  dirichlet_tags = DIRICHLET_tags, conformity = :H1)
Q = TestFESpace(Ω, reffeₚ,   conformity= :C0)
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

# Linear Solver
to_logfile("Stokes solve")
# u_ST, p_ST = solve(ls₀,stokes_op)
u_ST, p_ST = solve(stokes_op)

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
h2= get_cell_measure(Ω_f)
h = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω_f)),Ω_f)
τₘ = 1/(c₁*ν_f/h2 + c₂*(meas∘uₙₕ)/h)
τc = cc *(h2/(c₁*τₘ))

# Weak form
c(a,u,v) = 0.5*((∇(u)'⋅a)⋅v - u⋅(∇(v)'⋅a))
res(t,(u,p),(v,q)) = ∫( ∂t(u)⋅v  + c(u,u,v) + ε(v) ⊙ (σ_dev_f ∘ ε(u)) - p*(∇⋅v) + (∇⋅u)*q +
                        τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅u)) + τc*((∇⋅u)*(∇⋅v)) )dΩ_f +
                     ∫( 0.5*(u⋅v)*(u⋅n_Γout) )dΓout
jac(t,(u,p),(du,dp),(v,q)) = ∫( c(du,u,v) + c(u,du,v) + ε(v) ⊙ (σ_dev_f ∘ ε(du)) - dp*(∇⋅v) + (∇⋅du)*q +
τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅du) + (∇(du)'⋅u + ∇(u)'⋅du)⋅(∇(v)'⋅u)) +
τc*((∇⋅du)*(∇⋅v)) )dΩ_f +
∫( 0.5*((du⋅v)*(u⋅n_Γout)+(u⋅v)*(du⋅n_Γout)) )dΓout
jac_t(t,(u,p),(dut,dpt),(v,q)) = ∫( dut⋅v )dΩ_f

# Orthogonal projection
aη(η,κ) = ∫( τₘ*(η⋅κ) )dΩ_f
bη(κ) = ∫( τₘ*((∇(uₙₕ)'⋅uₙₕ)⋅κ) )dΩ_f
op_proj = AffineFEOperator(aη,bη,U(0.0),V)
# ls_proj = PETScLinearSolver()

# NS operator
op = TransientFEOperator(res,jac,jac_t, X, Y)
# op = TransientFEOperator(res,X, Y)

# Nonlinear Solver
# nls = PETScNonlinearSolver()

# Nonlinear Solver
ls = LUSolver()
nls = NLSolver(ls,show_trace=true,method=:newton,iterations=10)

# ODE solver
Δt = 0.1
tf = 0.1
ode_solver = GeneralizedAlpha(nls,Δt,ρ∞)

xₜ = solve(ode_solver,op,(xh₀,vh₀),t₀,tf)

# Postprocess
global tout = 0
createpvd("NS_test") do pvd
  for ((uh,ph),t) in xₜ
    to_logfile("Time: $t")
    to_logfile("=======================")
    Fx, Fy = sum(∫((n_ΓS ⋅ σ_dev_f(ε(uh))) - ph * n_ΓS) * dΓₛ)
    to_forcefile(t,Fx,Fy)
    if t>tout
      pvd[t] = createvtk(Ω,"NS_test_$t",cellfields=["u"=>uh,"p"=>ph,"un"=>uₙₕ,"eta_n"=>ηₙₕ])
      tout=t+Δtout
    end
    uₙₕ = interpolate!(uh,fv_u,U(t))
    ηₙₕ = solve(op_proj)
  end
end

close(io)
close(io_force)

end
