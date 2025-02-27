const g = 9.8
const treshold = 1e-6 # treshold for dealing with dry water states (h=0)

# DEFINITION OF SAINT VENANT FUNCTIONS 

struct SaintVenant <: AbstractEquationFun end

function flux(::SaintVenant, v)
    res = similar(v)
    g_half = g * 0.5 
    h = view(v, :, 1)
    for i in eachindex(h)
        if h[i] > treshold
            hu = v[i, 2]
            res[i, 1] = hu
            res[i, 2] = hu^2 / h[i] + g_half * h[i]^2
        else
            res[i, 1] = zero(eltype(v))
            res[i, 2] = zero(eltype(v))
        end
    end
    res
end

function Dflux(::SaintVenant, u)
end

# COMPUTING CFL CONDITION FOR SAINT VENANT

mutable struct CFLCacheSaintVenant <: CFLCacheType
    cfl::Float64
    eigenmax::Vector{Float64}
    function CFLCacheSaintVenant(uinit)
        eigenmax = zero(uinit[:,1])
        h = view(uinit, :, 1)
    for i in eachindex(h)
        eigenmax[i] = h[i] > treshold ? abs(uinit[i,2] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
    end
        new(zero(Float64), eigenmax)
    end
end

init_cfl_cache(::OneD, ::System, ::SaintVenant, equation, uinit) = CFLCacheSaintVenant(uinit)
function update_cflcache!(::OneD, ::System, ::SaintVenant, integrator::Integrator)
    @unpack u, cache = integrator
    @unpack eigenmax = cache.cfl_cache
    h = view(u, :, 1)
    for i in eachindex(h)
        eigenmax[i] = h[i] > treshold ? abs(u[i,2] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
    end
end

function CFL_cond(::SaintVenant, v)
    max = 0.0
    lamb = 0.0
    h = view(v, :, 1)
    for i in eachindex(h)
        hi = h[i]
        hu = v[i,2]
        if hi > treshold
            lamb = abs(hu / hi) + sqrt(g * abs(hi))
            if lamb > max
                max = lamb
            end
        end
    end
    return max
end

function CFL_cond!(::SaintVenant, integrator::Integrator)

    @unpack cfl_cache = integrator.cache

    # @unpack uprev = integrator
    # integrator.cfl = 0.0
    # lamb = 0.0

    # for i in eachindex(view(uprev, :, 1))
    #     h = uprev[i,1]
    #     hu = uprev[i,2]
    #     if h > treshold
    #         lamb = abs(hu / h) + sqrt(g * abs(h))
    #         if lamb > integrator.cfl
    #             integrator.cfl = lamb
    #         end
    #     end
    # end
    cfl_cache.cfl = maximum(cfl_cache.eigenmax)
end

function CFL_local!(::OneD, ::SaintVenant, integrator::Integrator, j::Int)

    @unpack uprev, cache, space_cache = integrator
    # @unpack stencil = cache
    @unpack eigenmax = cache.cfl_cache
    @unpack Nx = integrator.params.mesh

    # space_cache.cfl_loc = maximum(view(eigenmax, stencil))
    space_cache.cfl_loc = max(eigenmax[j], eigenmax[mod1(j+1,Nx)])

    # space_cache.cfl_loc = 0.0
    # lamb = 0.0

    # for i in eachindex(stencil)
    #     h = uprev[i,1]
    #     hu = uprev[i,2]
    #     if h > treshold
    #         lamb = abs(hu / h) + sqrt(g * abs(h))
    #         if lamb > space_cache.cfl_loc
    #             space_cache.cfl_loc = lamb
    #         end
    #     end
    # end
end

# SOME INITIAL CONDITIONS

# function init_lake_at_rest(x::T, znum::T; c=one(eltype(x))) where T<:AbstractVector
#     # znum = z.(x)
#     v = zeros(eltype(x), (length(x), 2))
#     for i in eachindex(x)
#         v[i,1] = max(zero(eltype(x)), c - znum[i])
#         v[i,2] = zero(eltype(x))
#     end
#     return v
# end

function init_lake_at_rest(x::T, znum::T; c=one(eltype(x))) where T<:AbstractVector
    nvar = ndims(znum)+1
    v = zeros(eltype(x), (size(znum)..., nvar))
    indices = indices = [Colon() for i in 1:nvar-1]
    for r in 1:nvar-1
        vr = view(v, indices..., r)
        for i in eachindex(x)
            vr[i] = max(zero(eltype(x)), c - znum[i])
        end
    end
    vend = view(v, indices..., nvar)
    vend .= zero(eltype(x))
    return v
end

# SOURCE TERM FOR SAINT VENANT EQUATION

struct Pointwise <: SourceDiscretize end
struct HRDisc <: SourceDiscretize end

struct TopoSource{sourceDiscretizeType<:SourceDiscretize, F1<:Base.Callable, F2<:Base.Callable} <: AbstractSource
    z::F1
    Dz::F2
    source_discretize::sourceDiscretizeType
end

z(ts::TopoSource, x::AbstractVector) = ts.z.(x)
Dz(ts::TopoSource, x::AbstractVector) = ts.Dz.(x)
Dz(ts::TopoSource, x::Real) = ts.Dz(x)

z(ts::TopoSource, mesh::OneDMesh) = ts.z.(mesh.x)
Dz(ts::TopoSource, mesh::OneDMesh) = ts.Dz.(mesh.x)

struct TopoSourceCache{znumType, DznumType} <: sourceCacheType
    znum::znumType
    Dznum::DznumType
    function TopoSourceCache(topo_source::TopoSource, source_discretize::SourceDiscretize, mesh::Mesh)
        znum = z(topo_source, mesh)
        Dznum = init_Dznum(source_discretize, topo_source, mesh)
        new{typeof(znum), typeof(Dznum)}(znum, Dznum)
    end
end

# INIT CACHE CONTENT
init_Dznum(::Pointwise, topo_source::TopoSource, mesh::Mesh) = Dz(topo_source, mesh)
init_Dznum(::HRDisc, args...) = nothing

# BY DEFAULT, POINTWISE DISCRETIZATION

# function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, mesh::OneDMesh)
#     @unpack x = mesh
#     res = zero(v)
#     for i in eachindex(view(v, :, 1))
#         res[i,2] = - v[i,1] * g * Dz(topo_source, x[i])
#     end
#     res
# end

# function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, mesh::OneDMesh, source_cache::TopoSourceCache)
#     @unpack x = mesh
#     @unpack Dznum = source_cache
#     res = similar(v)
#     for i in eachindex(view(v, :, 1))
#         res[i,1] = zero(eltype(v))
#         res[i,2] = - v[i,1] * g * Dznum[i]
#     end
#     res
# end

# function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, mesh::Mesh, source_cache::TopoSourceCache)
#     @unpack Dznum = source_cache
#     s = similar(v)
#     nvar = ndims(Dznum)+1
#     indices = [Colon() for _ in 1:nvar-1]
#     # FIRST EQUATION HAS ZERO SOURCE TERM
#     s1 = view(s, indices..., 1)
#     for i in eachindex(s1)
#         s1[i] = zero(eltype(v))
#     end
#     # OTHER EQUATIONS HAVE SPACE DERIVATED SOURCE TERMS
#     for r in 1:nvar-1
#         sr = view(v, indices..., r)
#         for i in eachindex(sr)
#             if nvar==2
#                 sr[i] = -v[i][1]*g*Dznum[i]
#             else
#                 sr[i] = -v[i][1]*g*Dznum[i][r]
#             end
#         end
#     end
#     s
# end

function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, mesh::Mesh, source_cache::TopoSourceCache)
    @unpack Dznum = source_cache
    s = similar(v)
    nvar = ndims(Dznum)+1
    # FIRST EQUATION HAS ZERO SOURCE TERM
    s1 = selectdim(s, nvar, 1)
    fill!(s1, zero(eltype(v)))
    # OTHER EQUATIONS HAVE SPACE DERIVATED SOURCE TERMS
    if nvar==2
        for i in eachindex(Dznum)
            s[i,2] = -v[i][1]*g*Dznum[i]
        end
    else
        for r in 2:nvar
            sr = selectdim(s, nvar, r)
            for i in eachindex(s1)
                sr[i] = -g * v[i, 1] * Dznum[i][r-1]
            end
        end
    end
    s
end

function discretize_sourceterm!(::Pointwise, integrator::Integrator)
    @unpack u, cache, source_cache = integrator
    @views h = u[:,1]
    for i in eachindex(h)
        cache.sourceterm[i,1] = zero(eltype(cache.sourceterm))
        cache.sourceterm[i,2] = -h[i]*g*source_cache.Dznum[i]
    end
end

# INITIALIZATION FUNCTIONS

function initialize_u(::OneD, ::EquationType, source::TopoSource, equation::AbstractEquation, params::Parameters)
    @unpack x = params.mesh
    znum = z(source, x)
    #(equation.initcond(x, znum), init_cache(source, source.source_discretize, x, znum))
    equation.initcond(x, znum)
end
init_sourceterm(source::TopoSource, args...) = discretize_sourceterm(source.source_discretize, source, args...)

# EXAMPLE OF CONFIGURATIONS FOR SAINT VENANT EQUATION

# 1 # LAKE AT REST WITH FLAT TOPOGRAPHY

zflat(x) = zero(x)
Dzflat(x) = zero(x)
FlatTopo = TopoSource(zflat, Dzflat, Pointwise())

# SaintVenantFlat = Equation(2, System(), SaintVenant(), x -> init_lake_at_rest(x,zflat), FlatTopo)
SaintVenantFlat = Equation(OneD(), 2, System(), SaintVenant(), (x,znum) -> init_lake_at_rest(x,znum), FlatTopo)

# 2 # LAKE AT REST WITH SINUSOIDAL TOPOGRAPHY

const freq = 1.0
const height = 0.5
zsinus(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
Dzsinus(x) = pi*freq*(sin(2*pi*freq * x))*height
BumpTopo = TopoSource(zsinus, Dzsinus, Pointwise())

SaintVenantAtRest = Equation(OneD(), 2, System(), SaintVenant(), (x,znum) -> init_lake_at_rest(x,znum), BumpTopo)