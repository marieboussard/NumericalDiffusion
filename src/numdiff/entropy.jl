abstract type AbstractEntropyFun end

struct EntropyFun{F1, F2} <: AbstractEntropyFun
    eta!::F1
    G!::F2
end

eta!(entfun::EntropyFun, args...) = entfun.eta!(args...)
G!(entfun::EntropyFun, args...) = entfun.eta!(args...)

# SOME CLASSICAL ENTROPIES

struct BurgersEnt <: AbstractEntropyFun end

eta!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = u^2
G!(::BurgersEnt, u::AbstractArray, res::AbstractArray) = @. res = 2*u^3/3

entropy(::Type{Burgers}) = BurgersEnt()