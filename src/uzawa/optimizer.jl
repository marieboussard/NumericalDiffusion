struct OptimizerOptions
    maxiter::Int
    eps::Float64
end

mutable struct OptimizerCache{ptype<:AbstractVector, gtype<:AbstractVector}
    Atp::ptype
    WAtp::ptype
    Agamma::gtype
    function OptimizerCache(p0::AbstractVector, Gc::AbstractVector)
        new{typeof(p0), typeof(Gc)}(zero(p0), zero(p0), zero(Gc))
    end
end


mutable struct Optimizer{wtype<:AbstractMatrix, atype<:AbstractMatrix, ptype<:AbstractVector, gtype<:AbstractVector, cachetype<:OptimizerCache}

    # Problem components
    W::wtype # Weights matrix
    A::atype # Constraints matrix
    b::ptype # Constraint vector 
    Gc::gtype

    # Variables
    p0::ptype
    gammaprev::gtype
    pprev::ptype
    gamma::gtype
    p::ptype

    niter::Int
    iterate_gap::Float64
    opts::OptimizerOptions
    mu::Float64
    cache::cachetype

    function Optimizer(Gc::AbstractVector, A::AbstractMatrix, b::AbstractVector; p0::AbstractVector=zero(b), W::AbstractMatrix=Matrix{eltype(Gc)}(I,length(Gc),length(Gc)), maxiter::Int=100, eps::Float64=1e-3)

        # Variables
        gammaprev = zero(Gc)
        pprev = zero(p0)
        gamma = zero(Gc)
        p = zero(p0)

        niter = 0
        iterate_gap = Inf
        opts = OptimizerOptions(maxiter, eps)
        mu = zero(Float64)
        cache = OptimizerCache(p0, Gc)

        new{typeof(W), typeof(A), typeof(p0), typeof(Gc), typeof(cache)}(W, A, b, Gc, p0, gammaprev, pprev, gamma, p, niter, iterate_gap, opts, mu, cache)
    end
end