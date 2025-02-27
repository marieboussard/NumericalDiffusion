struct SaintVenant2D <: AbstractEquationFun end

function flux_f(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    h = view(v, :, :, 1)
    hu = view(v, :, :, 2)
    hv = view(v, :, :, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            res[I, 1] = hu[I]
            res[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            res[I, 3] = hu[I]*hv[I] / h[I]
        else
            res[I, 1] = zero(eltype(v))
            res[I, 2] = zero(eltype(v))
            res[I, 3] = zero(eltype(v))
        end
    end
    res
end

function flux_f!(::SaintVenant2D, integrator)
    @unpack fcont = integrator.fcont
    @unpack u = integrator
    g_half = g * 0.5 
    h = view(u, :, :, 1)
    hu = view(u, :, :, 2)
    hv = view(u, :, :, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            fcont[I, 1] = hu[I]
            fcont[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
            fcont[I, 3] = hu[I]*hv[I] / h[I]
        else
            fcont[I, 1] = zero(eltype(u))
            fcont[I, 2] = zero(eltype(u))
            fcont[I, 3] = zero(eltype(u))
        end
    end
end

function flux_h(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    h = view(v, :, :, 1)
    hu = view(v, :, :, 2)
    hv = view(v, :, :, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            # hu = v[i, 2]
            # hv = v[i, 3]
            res[I, 1] = hv[I]
            res[I, 2] = hu[I]*hv[I] / h[I]
            res[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            res[I, 1] = zero(eltype(v))
            res[I, 2] = zero(eltype(v))
            res[I, 3] = zero(eltype(v))
        end
    end
    res
end

function flux_h!(::SaintVenant2D, integrator)
    @unpack hcont = integrator.fcont
    @unpack u = integrator
    g_half = g * 0.5 
    h = view(u, :, :, 1)
    hu = view(u, :, :, 2)
    hv = view(u, :, :, 3)
    for I in CartesianIndices(h)
        if h[I] > treshold
            hcont[I, 1] = hv[I]
            hcont[I, 2] = hu[I]*hv[I] / h[I]
            hcont[I, 3] = hv[I]^2 / h[I] + g_half * h[I]^2
        else
            hcont[I, 1] = zero(eltype(u))
            hcont[I, 2] = zero(eltype(u))
            hcont[I, 3] = zero(eltype(u))
        end
    end
end

# COMPUTING CFL CONDITION FOR SAINT VENANT 2D

mutable struct CFLCacheSaintVenant2D <: CFLCacheType
    cflx::Float64
    cfly::Float64
    xeigenmax::Matrix{Float64}
    yeigenmax::Matrix{Float64}
    function CFLCacheSaintVenant2D(uinit)
        xeigenmax = zero(uinit[:,:,1])
        yeigenmax = zero(uinit[:,:,1])
        h = view(uinit, :, :, 1)
        hu = view(uinit, :, :, 2)
        hv = view(uinit, :, :, 3)
    for i in eachindex(h)
        xeigenmax[i] = h[i] > treshold ? abs(hu[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
        yeigenmax[i] = h[i] > treshold ? abs(hv[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
    end
        new(zero(Float64), zero(Float64), xeigenmax, yeigenmax)
    end
end

init_cfl_cache(::TwoD, ::System, ::SaintVenant2D, equation, uinit) = CFLCacheSaintVenant2D(uinit)
function update_cflcache!(::TwoD, ::System, ::SaintVenant2D, integrator::Integrator)
    @unpack u, cache = integrator
    @unpack xeigenmax, yeigenmax = cache.cfl_cache
    h = view(u, :, :, 1)
    hu = view(u, :, :, 2)
    hv = view(u, :, :, 3)
    for i in eachindex(h)
        xeigenmax[i] = h[i] > treshold ? abs(hu[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(xeigenmax))
        yeigenmax[i] = h[i] > treshold ? abs(hv[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(yeigenmax))
    end
end

function CFL_cond2D!(::SaintVenant2D, integrator::Integrator)
    @unpack cfl_cache = integrator.cache
    cfl_cache.cflx = maximum(cfl_cache.xeigenmax)
    cfl_cache.cfly = maximum(cfl_cache.yeigenmax)
end

function CFL_local!(::TwoD, ::SaintVenant2D, integrator::Integrator, j::Int, k::Int)
    @unpack uprev, cache, space_cache = integrator
    @unpack xeigenmax, yeigenmax = cache.cfl_cache
    @unpack Nx, Ny = integrator.params.mesh
    space_cache.cflx_loc = max(xeigenmax[j, k], xeigenmax[mod1(j+1,Nx)])
    space_cache.cfly_loc = max(yeigenmax[j, k], yeigenmax[j, mod1(k+1,Ny)])
end

# SOME INITIAL CONDITIONS

# 1 # NO SOURCE TERM

function init_sv(x, y, args...)
    # return ones(eltype(x), (length(x), length(y), 2))
    return [1.0, 0.0, 0.0]
end

SaintVenant2Flat = Equation(TwoD(), 3, System(), SaintVenant2D(), init_sv)

# SOURCE TERM

z(ts::TopoSource, x::AbstractVector, y::AbstractVector) = [ts.z(xj, yk) for xj in x for yk in y]
Dz(ts::TopoSource, x::AbstractVector, y::AbstractVector) = [ts.Dz(xj, yk) for xj in x for yk in y]
z(ts::TopoSource, x::Real, y::Real) = ts.z(x, y)
Dz(ts::TopoSource, x::Real, y::Real) = ts.Dz(x, y)

function z(ts::TopoSource, mesh::TwoDMesh)
    @unpack x, y, Nx, Ny = mesh
    res = zeros(eltype(x), Nx, Ny)
    for j in 1:Nx
        for k in 1:Ny
            res[j,k] = ts.z(x[j], y[k])
        end
    end
end
function Dz(ts::TopoSource, mesh::TwoDMesh)
    @unpack x, y, Nx, Ny = mesh
    res = zeros(eltype(x), Nx, Ny, 2)
    for j in 1:Nx
        for k in 1:Ny
            dz = ts.Dz(x[j], y[k])
            res[j,k,1] = dz[1]
            res[j,k,2] = dz[2]
        end
    end
end

# BY DEFAULT, POINTWISE DISCRETIZATION

# function discretize_sourceterm(::Pointwise, topo_source::TopoSource , v, mesh::Mesh, source_cache::TopoSourceCache)
#     @unpack Dznum = source_cache
#     s = similar(v)
#     nvar = ndims(Dznum)+1
#     indices = indices = [Colon() for _ in 1:nvar-1]
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

# INITIALIZATION FUNCTIONS

# function initialize_u(::OneD, ::EquationType, source::TopoSource, equation::AbstractEquation, params::Parameters)
#     @unpack x = params.mesh
#     znum = z(source, x)
#     #(equation.initcond(x, znum), init_cache(source, source.source_discretize, x, znum))
#     equation.initcond(x, znum)
# end