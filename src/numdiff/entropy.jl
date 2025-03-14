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

struct BurgersEnt <: AbstractEntropyFun end

G(::BurgersEnt, u::Float64) = 2.0*u^3/3.0

# eta!(::BurgersEnt, u::Float64, res::Float64) = res = u^2
eta!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = u^2
# G!(::BurgersEnt, u::Float64, res::Float64) = res = 2*u^3/3
G!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = 2*u^3/3


entropy(::Type{Burgers}) = BurgersEnt()