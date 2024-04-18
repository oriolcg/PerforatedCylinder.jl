function get_stabilization_parameters(Ω,ν,c₁,c₂,cc)

  h = CellField(get_cell_measure(Ω),Ω)
  h2 = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω)),Ω)
  abs_(u) = (u⋅u)^(1/2)+1.0e-14
  dabs_(u,du) = (u⋅du)/abs_(u)
  τₘ⁻¹(u) = (c₁*ν/h2 + c₂*(abs_∘u)/h)
  τₘ(u) = 1/τₘ⁻¹(u)
  τc(u) = cc *(h2/(c₁*τₘ(u)))
  dτₘ(u,du) = -1.0/(τₘ⁻¹(u)*τₘ⁻¹(u)) * (c₂*(dabs_∘(u,du)))
  dτc(u,du) = -cc*h2/c₁ * (1/(τₘ(u)*τₘ(u))) * dτₘ(u,du)

  return τₘ, τc, dτₘ, dτc

end

# Stabilization operators
conv(a,∇u) = (∇u'⋅a)
ℒ(a,u,p) = (conv∘(a,∇(u))) + ∇(p)
∂ₐℒ(da,u) = conv∘(da,∇(u))
𝒫(a,u,p,η) = ℒ(a,u,p)-η
∂ₐ𝒫(da,u) = ∂ₐℒ(da,u)
uₛ(τₘ,a,u,p,η) = τₘ(a)*𝒫(a,u,p,η)
# ∂uₛ(a,u,p,η,da,du,dp,dη) = dτₘ(a,da)*𝒫(a,u,p,η) + τₘ(a)*(𝒫(a,du,dp,dη)+∂ₐ𝒫(da,u))

# Auxiliar functions
neg(a) = min(a,0.0)
conv_bc(a,∇u) = lazy_map(BroadcastingFieldOpMap(⋅), a, ∇u)


# Residual functions
# ==================

function f_conv(a,u,v,x)
  ax = (0.5*a)(x)
  ux = u(x)
  vx = v(x)
  ∇ux = ∇(u)(x)
  ∇vx = ∇(v)(x)
  cux = conv_bc(ax,∇ux)
  cvx = conv_bc(ax,∇vx)
  aux1 = lazy_map(BroadcastingFieldOpMap(⋅),cux,vx)
  aux2 = lazy_map(BroadcastingFieldOpMap(⋅),ux,cvx)
  return lazy_map(BroadcastingFieldOpMap(-),aux1,aux2)
end

function f_stab(τₘ,a,u,p,η,v,q,κ,x)
  # Functions
  ax = a(x)
  ηx = η(x)
  κx = κ(x)

  # Gradients
  ∇ux = ∇(u)(x)
  ∇px = ∇(p)(x)
  ∇vx = ∇(v)(x)
  ∇qx = ∇(q)(x)

  # Stabilization parameters
  τₘx = τₘ(a)(x)

  # Convective Terms
  cᵤx = conv_bc(ax,∇ux)
  cᵥx = conv_bc(ax,∇vx)

  # Operators
  ℒᵤx = lazy_map(BroadcastingFieldOpMap(+),cᵤx,∇px)
  𝒫ᵤx = lazy_map(BroadcastingFieldOpMap(-),ℒᵤx,ηx)
  ℒᵥx = lazy_map(BroadcastingFieldOpMap(+),cᵥx,∇qx)
  𝒫ᵥx = lazy_map(BroadcastingFieldOpMap(-),ℒᵥx,κx)

  # Subscales
  uₛx = lazy_map(BroadcastingFieldOpMap(*),τₘx,𝒫ᵤx)

  return lazy_map(BroadcastingFieldOpMap(⋅),uₛx,𝒫ᵥx)

end

function f_graddiv(τc,a,u,v,x)
  τcx = τc(a)(x)
  divux = (∇⋅u)(x)
  divvx = (∇⋅v)(x)

  divudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divux,divvx)
  return lazy_map(BroadcastingFieldOpMap(*),τcx,divudivvx)

end

# residual terms wrappers
lap(ν,u,v,dΩ) = ∫( 2ν*(ε(v) ⊙ ε(u)) )dΩ
div(u,q,dΩ) = ∫( q*(∇⋅u) )dΩ
cΓ(a,u,v,nΓ,dΓ) = ∫( (a⋅v)*(0.5*(u⋅nΓ)-neg∘(u⋅nΓ)) )dΓ
conv(a,u,v,dΩ) = own_integrate(x->f_conv(a,u,v,x),dΩ.quad)
stab(τₘ,a,u,p,η,v,q,κ,dΩ) = own_integrate(x->f_stab(τₘ,a,u,p,η,v,q,κ,x),dΩ.quad)
graddiv(τc,a,u,v,dΩ) = own_integrate(x->f_graddiv(τc,a,u,v,x), dΩ.quad)

# Jacobian functions
# ==================
function f_dconv(a,u,da,du,v,x)
  ax = (0.5*a)(x)
  ux = u(x)
  vx = v(x)
  dax = (0.5*da)(x)
  dux = du(x)
  ∇ux = ∇(u)(x)
  ∇vx = ∇(v)(x)
  ∇dux = ∇(du)(x)

  cadux = conv_bc(ax,∇dux)
  cavx = conv_bc(ax,∇vx)
  cdaux = conv_bc(dax,∇ux)
  cdavx = conv_bc(dax,∇vx)

  aux1 = lazy_map(BroadcastingFieldOpMap(+),cadux,cdaux)
  aux2 = lazy_map(BroadcastingFieldOpMap(⋅),aux1,vx)
  aux3 = lazy_map(BroadcastingFieldOpMap(⋅),dux,cavx)
  aux4 = lazy_map(BroadcastingFieldOpMap(⋅),ux,cdavx)
  aux4 = lazy_map(BroadcastingFieldOpMap(+),aux3,aux4)

  return lazy_map(BroadcastingFieldOpMap(-),aux2,aux4)
end

function f_dstab(τₘ,dτₘ,a,u,p,η,da,du,dp,dη,v,q,κ,x)

  # Functions
  ax = a(x)
  ηx = η(x)
  κx = κ(x)
  dax = da(x)
  dηx = dη(x)

  # Gradients
  ∇ux = ∇(u)(x)
  ∇px = ∇(p)(x)
  ∇vx = ∇(v)(x)
  ∇qx = ∇(q)(x)
  ∇dux = ∇(du)(x)
  ∇dpx = ∇(dp)(x)

  # Stabilization parameters
  τₘx = τₘ(a)(x)
  dτₘx = dτₘ(a,da)(x)

  # Convective Terms
  cᵤx = conv_bc(ax,∇ux)
  cᵥx = conv_bc(ax,∇vx)
  cdux = conv_bc(ax,∇dux)

  # Operators
  ℒᵤx = lazy_map(BroadcastingFieldOpMap(+),cᵤx,∇px)
  𝒫ᵤx = lazy_map(BroadcastingFieldOpMap(-),ℒᵤx,ηx)
  ℒdux = lazy_map(BroadcastingFieldOpMap(+),cdux,∇dpx)
  𝒫dux = lazy_map(BroadcastingFieldOpMap(-),ℒdux,dηx)
  ∂ₐ𝒫x = conv_bc(dax,∇ux)
  d𝒫ᵤx = lazy_map(BroadcastingFieldOpMap(+),𝒫dux,∂ₐ𝒫x)
  ℒᵥx = lazy_map(BroadcastingFieldOpMap(+),cᵥx,∇qx)
  𝒫ᵥx = lazy_map(BroadcastingFieldOpMap(-),ℒᵥx,κx)

  # Subscales
  uₛx = lazy_map(BroadcastingFieldOpMap(*),τₘx,𝒫ᵤx)
  ∂uₛx = lazy_map(BroadcastingFieldOpMap(+),
    lazy_map(BroadcastingFieldOpMap(*),dτₘx,𝒫ᵤx),
    lazy_map(BroadcastingFieldOpMap(*),τₘx,d𝒫ᵤx))
  ∂ₐ𝒫ᵥx = conv_bc(dax,∇vx)

  aux1 = lazy_map(BroadcastingFieldOpMap(⋅),∂uₛx,𝒫ᵥx)
  aux2 = lazy_map(BroadcastingFieldOpMap(⋅),uₛx,∂ₐ𝒫ᵥx)

  return lazy_map(BroadcastingFieldOpMap(+),aux1,aux2)

end

function f_dgraddiv(τc,dτc,a,u,da,du,v,x)
  τcx = τc(a)(x)
  divux = (∇⋅u)(x)
  divvx = (∇⋅v)(x)
  dτcx = dτc(a,da)(x)
  divdux = (∇⋅du)(x)

  divdudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divdux,divvx)
  τdivx = lazy_map(BroadcastingFieldOpMap(*),τcx,divdudivvx)
  divudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divux,divvx)
  dτdivx = lazy_map(BroadcastingFieldOpMap(*),dτcx,divudivvx)

  return lazy_map(BroadcastingFieldOpMap(+),τdivx,dτdivx)

end

function f_div(du,dp,v,q,x)
  divux = (∇⋅du)(x)
  divvx = (∇⋅v)(x)
  dpx = dp(x)
  qx = q(x)
  qdivu = lazy_map(BroadcastingFieldOpMap(*),divux,qx)
  pdivv = lazy_map(BroadcastingFieldOpMap(*),divvx,dpx)
  return lazy_map(BroadcastingFieldOpMap(-),qdivu,pdivv)
end

# Jacobian terms wrappers
dconv(a,u,da,du,v,dΩ) = own_integrate(x->f_dconv(a,u,da,du,v,x),dΩ.quad)
dstab(τₘ,dτₘ,a,u,p,η,da,du,dp,dη,v,q,κ,dΩ) = own_integrate(x->f_dstab(τₘ,dτₘ,a,u,p,η,da,du,dp,dη,v,q,κ,x),dΩ.quad)
dgraddiv(τc,dτc,a,u,da,du,v,dΩ) = own_integrate(x->f_dgraddiv(τc,dτc,a,u,da,du,v,x), dΩ.quad)

# Integration function
# ====================
function own_integrate(f::Function,quad::CellQuadrature)
  x = get_cell_points(quad)
  bx = f(x)
  if quad.data_domain_style == PhysicalDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    result = lazy_map(IntegrationMap(),bx,quad.cell_weight)
  elseif quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == PhysicalDomain()
    cell_map = get_cell_map(quad.trian)
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    result = lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)
  elseif quad.data_domain_style == ReferenceDomain() &&
            quad.integration_domain_style == ReferenceDomain()
    cell_map = Fill(GenericField(identity),length(bx))
    cell_Jt = lazy_map(∇,cell_map)
    cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
    result = lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)
  else
    @notimplemented
  end
  c__ = DomainContribution()
  add_contribution!(c__,quad.trian,result)
  return c__
end
