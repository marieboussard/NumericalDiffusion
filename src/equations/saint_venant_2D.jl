struct SaintVenant2D <: AbstractEquationFun end

function flux_f(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    h = view(v, :, 1)
    for i in eachindex(h)
        if h[i] > treshold
            hu = v[i, 2]
            hv = v[i, 3]
            res[i, 1] = hu
            res[i, 2] = hu^2 / h[i] + g_half * h[i]^2
            res[i, 3] = hu*hv / h[i]
        else
            res[i, 1] = zero(eltype(v))
            res[i, 2] = zero(eltype(v))
            res[i, 3] = zero(eltype(v))
        end
    end
    res
end

function flux_h(::SaintVenant2D, v)
    res = similar(v)
    g_half = g * 0.5 
    h = view(v, :, 1)
    for i in eachindex(h)
        if h[i] > treshold
            hu = v[i, 2]
            hv = v[i, 3]
            res[i, 1] = hv
            res[i, 2] = hu*hv / h[i]
            res[i, 3] = hv^2 / h[i] + g_half * h[i]^2
        else
            res[i, 1] = zero(eltype(v))
            res[i, 2] = zero(eltype(v))
            res[i, 3] = zero(eltype(v))
        end
    end
    res
end

# COMPUTING CFL CONDITION FOR SAINT VENANT 2D

mutable struct CFLCacheSaintVenant2D <: CFLCacheType
    cflx::Float64
    cfly::Float64
    xeigenmax::Vector{Float64}
    yeigenmax::Vector{Float64}
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
    hu = view(uinit, :, :, 2)
    hv = view(uinit, :, :, 3)
    for i in eachindex(h)
        xeigenmax[i] = h[i] > treshold ? abs(hu[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
        yeigenmax[i] = h[i] > treshold ? abs(hv[i] / h[i]) + sqrt(g*h[i]) : zero(eltype(eigenmax))
    end
end

function CFL_cond2D!(::SaintVenant2D, integrator::Integrator)
    @unpack cfl_cache = integrator.cache
    cfl_cache.cflx = maximum(cfl_cache.xeigenmax)
    cfl_cache.cfly = maximum(cfl_cache.yeigenmax)
end

function CFL_local2D!(::SaintVenant2D, integrator::Integrator, j::Int)
    @unpack uprev, cache, space_cache = integrator
    @unpack xeigenmax, yeigenmax = cache.cfl_cache
    @unpack Nx, Ny = integrator.params.mesh
    space_cache.cflx_loc = max(xeigenmax[j], xeigenmax[mod1(j+1,Nx)])
    space_cache.cfly_loc = max(yeigenmax[j], yeigenmax[mod1(j+1,Ny)])
end

# SOME INITIAL CONDITIONS

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
    vend = view(v, indices, nvar)
    vend .= zero(eltype(x))
    return v
end