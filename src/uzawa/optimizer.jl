struct OptimizerOptions
    maxiter::Int
    eps::Float64
    start_with_gamma::Bool
end

mutable struct OptimizerCache{wtype<:AbstractMatrix, ptype<:AbstractVector, gtype<:AbstractVector}
    Winv::wtype
    Atp::gtype
    WAtp::gtype
    Agamma::ptype
    function OptimizerCache(W::AbstractMatrix, p0::AbstractVector, Gc::AbstractVector)
        Winv = zero(W)
        for j in 1:size(W)[1]
            Winv[j,j] = 1.0/W[j,j]
        end
        new{typeof(Winv), typeof(p0), typeof(Gc)}(Winv, zero(Gc), zero(Gc), zero(p0))
    end
end


mutable struct Optimizer{wtype<:AbstractMatrix, atype<:AbstractMatrix, ptype<:AbstractVector, gtype<:AbstractVector, cachetype<:OptimizerCache}

    # Problem components
    W::wtype # Weights matrix
    A::atype # Constraints matrix
    b::ptype # Constraint vector 
    Gc::gtype

    # Variables
    gamma0::gtype
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

    function Optimizer(Gc::AbstractVector, A::AbstractMatrix, b::AbstractVector; p0::AbstractVector=zero(b), gamma0=zero(Gc)::AbstractVector, W::AbstractMatrix=Matrix{eltype(Gc)}(I,length(Gc),length(Gc)), maxiter::Int=100, eps::Float64=1e-5, start_with_gamma::Bool=false)

        # Variables
        gammaprev = zero(Gc)
        pprev = copy(p0)
        gamma = zero(Gc)
        p = zero(p0)

        niter = 0
        iterate_gap = Inf
        opts = OptimizerOptions(maxiter, eps, start_with_gamma)
        mu = zero(Float64)
        cache = OptimizerCache(W, p0, Gc)

        new{typeof(W), typeof(A), typeof(p0), typeof(Gc), typeof(cache)}(W, A, b, Gc, gamma0, p0, gammaprev, pprev, gamma, p, niter, iterate_gap, opts, mu, cache)
    end
end