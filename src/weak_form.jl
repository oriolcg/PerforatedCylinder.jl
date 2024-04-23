function get_stabilization_parameters(ν,c₁,c₂,cc)

  abs_(u) = (u⋅u)^(1/2)+1.0e-14
  dabs_(u,du) = (u⋅du)/abs_(u)
  τₘ⁻¹(u,h,h2) = (c₁*ν/h2 + c₂*(abs_∘u)/h)
  τₘ(u,h,h2) = 1/τₘ⁻¹(u,h,h2)
  τc(u,h,h2) = cc *(h2/(c₁*τₘ(u,h,h2)))
  dτₘ(u,du,h,h2) = -1.0/(τₘ⁻¹(u,h,h2)*τₘ⁻¹(u,h,h2)) * (c₂*(dabs_∘(u,du)))
  dτc(u,du,h,h2) = -cc*h2/c₁ * (1/(τₘ(u,h,h2)*τₘ(u,h,h2))) * dτₘ(u,du,h,h2)

  return τₘ, τc, dτₘ, dτc

end
function get_stabilization_parameters_(ν,c₁,c₂,cc)

  abs_(u) = (u⋅u).^(1/2)+1.0e-14
  dabs_(u,du) = (u⋅du)/abs_(u)
  τₘ⁻¹(u,h,h2) = (c₁*ν/h2 + c₂*(abs_(u))/h)
  τₘ(u,h,h2) = 1/τₘ⁻¹(u,h,h2)
  τc(u,h,h2) = cc *(h2/(c₁*τₘ(u,h,h2)))
  dτₘ(u,du,h,h2) = -1.0/(τₘ⁻¹(u,h,h2)*τₘ⁻¹(u,h,h2)) * (c₂*(dabs_(u,du)))
  dτc(u,du,h,h2) = -cc*h2/c₁ * (1/(τₘ(u,h,h2)*τₘ(u,h,h2))) * dτₘ(u,du,h,h2)

  return τₘ, τc, dτₘ, dτc

end

function get_mesh_sizes(Ω)
  h = CellField(get_cell_measure(Ω),Ω)
  h2 = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω)),Ω)
  return h,h2
end

function get_mesh_sizes(Ω::GridapDistributed.DistributedTriangulation)
  h2map = map(Ω.trians) do trian
    CellField(get_cell_measure(trian),trian)
  end
  h2 = DistributedCellField(h2map,Ω)
  hmap = map(Ω.trians) do trian
    CellField(lazy_map(dx->dx^(1/2),get_cell_measure(trian)),trian)
  end
  h = DistributedCellField(hmap,Ω)
  return h,h2
end

# Stabilization operators
conv(a,∇u) = (∇u'⋅a)
ℒ(a,u,p) = (conv∘(a,∇(u))) + ∇(p)
∂ₐℒ(da,u) = conv∘(da,∇(u))
𝒫(a,u,p,η) = ℒ(a,u,p)-η
∂ₐ𝒫(da,u) = ∂ₐℒ(da,u)
uₛ(τₘ,h,h2,a,u,p,η) = τₘ(a,h,h2)*𝒫(a,u,p,η)
∂uₛ(τₘ,dτₘ,h,h2,a,u,p,η,da,du,dp,dη) = dτₘ(a,da,h,h2)*𝒫(a,u,p,η) + τₘ(a,h,h2)*(𝒫(a,du,dp,dη)+∂ₐ𝒫(da,u))

# Auxiliar functions
neg(a) = min(a,0.0)
conv_bc(a,∇u) = lazy_map(BroadcastingFieldOpMap(⋅), a, ∇u)
conv_(a,∇u) = (∇u'⋅a)


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

function f_stab(τₘ,h,h2,a,u,p,η,v,q,κ,x)
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
  τₘx = τₘ(a,h,h2)(x)

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

function f_stab_expl(τₘ,h,h2,a,u,p,η,v,q,x)
  # Functions
  ax = a(x)
  ηx = η(x)

  # Gradients
  ∇ux = ∇(u)(x)
  ∇px = ∇(p)(x)
  ∇vx = ∇(v)(x)
  ∇qx = ∇(q)(x)

  # Stabilization parameters
  τₘx = τₘ(a,h,h2)(x)

  # Convective Terms
  cᵤx = conv_bc(ax,∇ux)
  cᵥx = conv_bc(ax,∇vx)

  # Operators
  ℒᵤx = lazy_map(BroadcastingFieldOpMap(+),cᵤx,∇px)
  𝒫ᵤx = lazy_map(BroadcastingFieldOpMap(-),ℒᵤx,ηx)
  ℒᵥx = lazy_map(BroadcastingFieldOpMap(+),cᵥx,∇qx)

  # Subscales
  uₛx = lazy_map(BroadcastingFieldOpMap(*),τₘx,𝒫ᵤx)

  return lazy_map(BroadcastingFieldOpMap(⋅),uₛx,ℒᵥx)

end

function f_graddiv(τc,h,h2,a,u,v,x)
  τcx = τc(a,h,h2)(x)
  divux = (∇⋅u)(x)
  divvx = (∇⋅v)(x)

  divudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divux,divvx)
  return lazy_map(BroadcastingFieldOpMap(*),τcx,divudivvx)

end

# residual terms wrappers
lap(ν,u,v,dΩ) = ∫( 2ν*(ε(v) ⊙ ε(u)) )dΩ
div(u,q,dΩ) = ∫( q*(∇⋅u) )dΩ
cΓ(a,u,v,nΓ,dΓ) = ∫( (a⋅v)*(0.5*(u⋅nΓ)-neg∘(u⋅nΓ)) )dΓ
conv(a,u,v,dΩ::Measure) = own_integrate(x->f_conv(a,u,v,x),dΩ)
stab(τₘ,h,h2,a,u,p,η,v,q,κ,dΩ::Measure) = own_integrate(x->f_stab(τₘ,h,h2,a,u,p,η,v,q,κ,x),dΩ)
stab_expl(τₘ,h,h2,a,u,p,η,v,q,dΩ::Measure) = own_integrate(x->f_stab_expl(τₘ,h,h2,a,u,p,η,v,q,x),dΩ)
graddiv(τc,h,h2,a,u,v,dΩ::Measure) = own_integrate(x->f_graddiv(τc,h,h2,a,u,v,x), dΩ)
function conv(a,u,v,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(a.fields,u.fields,v.fields,dΩ.measures) do af,uf,vf,m
    conv(af,uf,vf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function stab(τₘ,h,h2,a,u,p,η,v,q,κ,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,p.fields,η.fields,v.fields,q.fields,κ.fields,dΩ.measures) do hf,h2f,af,uf,pf,ηf,vf,qf,κf,m
    stab(τₘ,hf,h2f,af,uf,pf,ηf,vf,qf,κf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function stab_expl(τₘ,h,h2,a,u,p,η,v,q,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,p.fields,η.fields,v.fields,q.fields,dΩ.measures) do hf,h2f,af,uf,pf,ηf,vf,qf,m
    stab_expl(τₘ,hf,h2f,af,uf,pf,ηf,vf,qf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function graddiv(τc,h,h2,a,u,v,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,v.fields,dΩ.measures) do hf,h2f,af,uf,vf,m
    graddiv(τc,hf,h2f,af,uf,vf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
# conv(a,u,v,dΩ::GridapDistributed.DistributedMeasure) = ∫(0.5*((conv_∘(a,∇(u)))⋅v - u⋅(conv_∘(a,∇(v)))))dΩ
# stab(τₘ,a,u,p,η,v,q,κ,dΩ::GridapDistributed.DistributedMeasure) = ∫( uₛ(τₘ,a,u,p,η)⋅𝒫(a,v,q,κ))dΩ
# graddiv(τc,a,u,v,dΩ::GridapDistributed.DistributedMeasure) = ∫( τc(a)*((∇⋅u)*(∇⋅v)) )dΩ

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
  aux5 = lazy_map(BroadcastingFieldOpMap(+),aux3,aux4)

  return lazy_map(BroadcastingFieldOpMap(-),aux2,aux5)
end

function f_dstab(τₘ,dτₘ,h,h2,a,u,p,η,da,du,dp,dη,v,q,κ,x)

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
  τₘx = τₘ(a,h,h2)(x)
  dτₘx = dτₘ(a,da,h,h2)(x)

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

function f_dstab_expl(τₘ,h,h2,a,u,p,η,du,dp,v,q,x)

  # Functions
  ax = a(x)
  ηx = η(x)

  # Gradients
  ∇ux = ∇(u)(x)
  ∇vx = ∇(v)(x)
  ∇qx = ∇(q)(x)
  ∇dux = ∇(du)(x)
  ∇dpx = ∇(dp)(x)

  # Stabilization parameters
  τₘx = τₘ(a,h,h2)(x)

  # Convective Terms
  cᵥx = conv_bc(ax,∇vx)
  cdux = conv_bc(ax,∇dux)

  # # Operators
  # ℒdux = lazy_map(BroadcastingFieldOpMap(+),cdux,∇px)
  # ℒᵥx = lazy_map(BroadcastingFieldOpMap(+),cᵥx,∇qx)

  # # Subscales
  # uₛx = lazy_map(BroadcastingFieldOpMap(*),τₘx,𝒫ᵤx)

  # return lazy_map(BroadcastingFieldOpMap(⋅),uₛx,ℒᵥx)

  # Operators
  ℒdux = lazy_map(BroadcastingFieldOpMap(+),cdux,∇dpx)
  # 𝒫dux = lazy_map(BroadcastingFieldOpMap(-),ℒdux,ηx)
  ℒᵥx = lazy_map(BroadcastingFieldOpMap(+),cᵥx,∇qx)

  # Subscales
  ∂uₛx = lazy_map(BroadcastingFieldOpMap(*),τₘx,ℒdux)

  return lazy_map(BroadcastingFieldOpMap(⋅),∂uₛx,ℒᵥx)

end

function f_dgraddiv(τc,dτc,h,h2,a,u,da,du,v,x)
  τcx = τc(a,h,h2)(x)
  divux = (∇⋅u)(x)
  divvx = (∇⋅v)(x)
  dτcx = dτc(a,da,h,h2)(x)
  divdux = (∇⋅du)(x)

  divdudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divdux,divvx)
  τdivx = lazy_map(BroadcastingFieldOpMap(*),τcx,divdudivvx)
  divudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divux,divvx)
  dτdivx = lazy_map(BroadcastingFieldOpMap(*),dτcx,divudivvx)

  return lazy_map(BroadcastingFieldOpMap(+),τdivx,dτdivx)

end

function f_dgraddiv_expl(τc,h,h2,a,u,du,v,x)
  τcx = τc(a,h,h2)(x)
  divvx = (∇⋅v)(x)
  divdux = (∇⋅du)(x)

  divdudivvx = lazy_map(BroadcastingFieldOpMap(⋅),divdux,divvx)

  return lazy_map(BroadcastingFieldOpMap(*),τcx,divdudivvx)

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
dconv(a,u,da,du,v,dΩ::Measure) = own_integrate(x->f_dconv(a,u,da,du,v,x),dΩ)
dstab(τₘ,dτₘ,h,h2,a,u,p,η,da,du,dp,dη,v,q,κ,dΩ::Measure) = own_integrate(x->f_dstab(τₘ,dτₘ,h,h2,a,u,p,η,da,du,dp,dη,v,q,κ,x),dΩ)
dstab_expl(τₘ,h,h2,a,u,p,η,du,dp,v,q,dΩ::Measure) = own_integrate(x->f_dstab_expl(τₘ,h,h2,a,u,p,η,du,dp,v,q,x),dΩ)
dgraddiv(τc,dτc,h,h2,a,u,da,du,v,dΩ::Measure) = own_integrate(x->f_dgraddiv(τc,dτc,h,h2,a,u,da,du,v,x), dΩ)
dgraddiv_expl(τc,h,h2,a,u,du,v,dΩ::Measure) = own_integrate(x->f_dgraddiv_expl(τc,h,h2,a,u,du,v,x), dΩ)
function dconv(a,u,da,du,v,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(a.fields,u.fields,da.fields,du.fields,v.fields,dΩ.measures) do af,uf,daf,duf,vf,m
    dconv(af,uf,daf,duf,vf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function dstab(τₘ,dτₘ,h,h2,a,u,p,η,da,du,dp,dη,v,q,κ,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,p.fields,η.fields,da.fields,du.fields,dp.fields,dη.fields,v.fields,q.fields,κ.fields,dΩ.measures) do hf,h2f,af,uf,pf,ηf,daf,duf,dpf,dηf,vf,qf,κf,m
    dstab(τₘ,dτₘ,hf,h2f,af,uf,pf,ηf,daf,duf,dpf,dηf,vf,qf,κf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function dstab_expl(τₘ,h,h2,a,u,p,η,du,dp,v,q,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,p.fields,η.fields,du.fields,dp.fields,v.fields,q.fields,dΩ.measures) do hf,h2f,af,uf,pf,ηf,duf,dpf,vf,qf,m
    dstab_expl(τₘ,hf,h2f,af,uf,pf,ηf,duf,dpf,vf,qf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function dgraddiv(τc,dτc,h,h2,a,u,da,du,v,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,da.fields,du.fields,v.fields,dΩ.measures) do hf,h2f,af,uf,daf,duf,vf,m
    dgraddiv(τc,dτc,hf,h2f,af,uf,daf,duf,vf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
function dgraddiv_expl(τc,h,h2,a,u,du,v,dΩ::GridapDistributed.DistributedMeasure)
  contribs = map(h.fields,h2.fields,a.fields,u.fields,du.fields,v.fields,dΩ.measures) do hf,h2f,af,uf,duf,vf,m
    dgraddiv_expl(τc,hf,h2f,af,uf,duf,vf,m)
  end
  GridapDistributed.DistributedDomainContribution(contribs)
end
# dconv(a,u,da,du,v,dΩ::GridapDistributed.DistributedMeasure) = conv(da,u,v,dΩ) + conv(a,du,v,dΩ)
# dstab(τₘ,dτₘ,a,u,p,η,da,du,dp,dη,v,q,κ,dΩ::GridapDistributed.DistributedMeasure) =
#     ∫( ∂uₛ(τₘ,dτₘ,a,u,p,η,da,du,dp,dη)⋅𝒫(a,v,q,κ) )dΩ +
#     ∫( uₛ(τₘ,a,u,p,η)⋅∂ₐ𝒫(da,v) )dΩ
# dgraddiv(τc,dτc,a,u,da,du,v,dΩ::GridapDistributed.DistributedMeasure) =
#     ∫( τc(a)*((∇⋅du)*(∇⋅v)) )dΩ +
#     ∫( dτc(a,da)*((∇⋅u)*(∇⋅v)) )dΩ

# Integration function
# ====================
function own_integrate(f::Function,dΩ::Measure)
  quad = dΩ.quad
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

# function own_integrate(f,dΩ::GridapDistributed.DistributedMeasure)
#   contribs = map(dΩ.measures) do m
#     own_integrate(f,m)
#   end
#   DistributedDomainContribution(contribs)
# end
