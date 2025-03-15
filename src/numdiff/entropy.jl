abstract type AbstractEntropyFun end

struct EntropyFun{F1, F2, F3} <: AbstractEntropyFun
    eta!::F1
    G!::F2
    G::F3
end

eta!(entfun::EntropyFun, args...) = entfun.eta!(args...)
G!(entfun::EntropyFun, args...) = entfun.G!(args...)
G(entfun::EntropyFun, args...) = entfun.G(args...)

# SOME CLASSICAL ENTROPIES

# BURGERS 
struct BurgersEnt <: AbstractEntropyFun end

G(::BurgersEnt, u::Float64) = 2.0*u^3/3.0
G(::BurgersEnt, u::AbstractVector) = 2.0*u[1]^3/3.0

# eta!(::BurgersEnt, u::Float64, res::Float64) = res = u^2
eta!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = u^2
# G!(::BurgersEnt, u::Float64, res::Float64) = res = 2*u^3/3
G!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = 2*u^3/3

entropy(::Type{Burgers}) = BurgersEnt()

# SAINT VENANT
struct SaintVenantEnt <: AbstractEntropyFun end
function eta!(::SaintVenantEnt, v::AbstractArray, res::AbstractArray)
    g_half = g * 0.5 
    # if ndims(v)==1
    #     if v[1] > treshold
    #         res[1] = v[2]
    #         res[2] = v[2]^2/v[1] + g_half * v[1]^2
    #     else
    #         fill!(res, zero(eltype(v)))
    #     end
    # else
    h = view(v, :, 1)
    for i in eachindex(h)
        if h[i] > treshold
            hu = v[i, 2]
            res[i] = 0.5*hu^2 / h[i] + g_half * h[i]^2
        else
            res[i] = zero(eltype(v))
        end
    end
    #end
end
function G!(::SaintVenantEnt, v::AbstractArray, res::AbstractArray)
    g_half = g * 0.5 
    # if ndims(v)==1
    #     if v[1] > treshold
    #         res[1] = v[2]
    #         res[2] = v[2]^2/v[1] + g_half * v[1]^2
    #     else
    #         fill!(res, zero(eltype(v)))
    #     end
    # else
    h = view(v, :, 1)
    for i in eachindex(h)
        if h[i] > treshold
            hu = v[i, 2]
            res[i] = (0.5*hu^2 / h[i] + g_half * h[i]^2)*hu/h[i]
        else
            res[i] = zero(eltype(v))
        end
    end
    #end
end
G(::SaintVenantEnt, v::AbstractVector) = v[1] > treshold ? (0.5*v[2]^2/v[1] + 0.5*g*v[1]^2)*v[2]/v[1] : zero(eltype(v))

# etasource!(::SaintVenantEnt, v::AbstractArray, sourceterm::AbstractArray, res::AbstractArray) = @. res += view(v,:,1)*g*sourceterm
etasource!(::SaintVenantEnt, v::AbstractArray, sourceterm::AbstractArray, res::AbstractArray) =  res .+= view(v,:,1).*g.*sourceterm
etasource!(entfun::SaintVenantEnt, v::AbstractArray, source_cache::SourceCache, res::AbstractArray) = etasource!(entfun, v, source_cache.znum, res)
G!(::SaintVenantEnt, v::AbstractArray, z::AbstractArray, res::AbstractArray) = @. res += view(v,:,1)*g*z
Gsource(::SaintVenantEnt, v::AbstractVector, z::Float64) = v[1]*g*z
entropy(::Type{SaintVenant}) = SaintVenantEnt()