mutable struct EstimatorCache{utype <: AbstractArray, cflCacheType<:CFLCache, ktype} <: Cache
    sL::Int
    sR::Int
    indices::Vector{Int}
    utilde::utype
    uhat::utype
    fcont_tilde::utype
    ftilde::utype
    cfl_cache::cflCacheType
    eta_tilde::utype
    eta_hat::utype
    K::ktype
    GK::Float64
    
    function EstimatorCache(equation::Equation, time_scheme::TimeScheme, space_scheme::SpaceScheme, u::AbstractArray)
        sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)
        indices = zeros(Int64, 2*(sL+sR))
        utilde = init_utilde(equation.dim, equation.eqtype, u, sL, sR)
        uhat = init_uhat(equation.dim, equation.eqtype, u, sL, sR)
        fcont_tilde = zero(utilde)
        ftilde = init_ftilde(equation.dim, equation.eqtype, u, sL, sR)
        cfl_cache = init_cfl_cache(equation.dim, equation.eqtype, equation.funcs, equation, utilde)
        eta_tilde = zero(uhat)
        eta_hat = zero(uhat)
        K = init_K(equation.dim, equation.eqtype)
        GK = zero(Float64)
        new{typeof(utilde), typeof(cfl_cache), typeof(K)}(sL, sR, indices, utilde, uhat, fcont_tilde, ftilde, cfl_cache, eta_tilde, eta_hat, K, GK)
    end
end

init_utilde(::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR))
init_uhat(::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR))
init_ftilde(::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)+1)
init_K(::OneD, ::Scalar) = zero(Float64)


mutable struct Estimator{methodType<:QuantifMethod, ecacheType<:EstimatorCache, mcacheType<:MethodCache, scacheType<:SpaceCache, tcacheType<:TimeCache, entfunType<:AbstractEntropyFun, etaType <: AbstractArray, diffType<:AbstractArray}

    # PROBLEM COMPONENTS
    sol::Solution

    # QUANTIFICATION METHOD
    method::methodType

    # CACHE
    cache::ecacheType
    method_cache::mcacheType
    space_cache::scacheType
    time_cache::tcacheType
    
    # ENTROPY
    entfun::entfunType
    etacont_init::etaType
    etacont::etaType

    # BOUNDS FOR NUMERICAL ENTROPY FLUX
    # Gcont::gType
    m::etaType
    M::etaType

    # RESULTS
    D::diffType

    function Estimator(sol::Solution, method::QuantifMethod; kwargs...)

        # INIT CACHE
        cache = EstimatorCache(sol.equation, sol.time_scheme, sol.space_scheme, sol.u)
        method_cache = init_cache(method, sol.equation, sol.u)
        space_cache = init_cache(sol.space_scheme)
        time_cache = init_cache(sol.time_scheme)
        
        # INIT ENTROPY
        entfun = entropy(typeof(sol.equation.funcs))
        etacont_init = init_etacont(sol.equation.dim, sol.equation.eqtype, sol.u)
        etacont = zero(etacont_init)

        # INIT NUMERICAL ENTROPY FLUX BOUNDS
        m = zero(etacont_init)
        M = zero(etacont_init)

        # INIT DIFFUSION
        D = zero(sol.uinit)

        new{typeof(method), typeof(cache), typeof(method_cache), typeof(space_cache), typeof(time_cache), typeof(entfun), typeof(etacont), typeof(D)}(sol, method, cache, method_cache, space_cache, time_cache, entfun, etacont_init, etacont, m, M, D)
    end
end

# FASTER ACCESS TO SOME FIELDS OF ESTIMATOR

function Base.getproperty(estimator::Estimator, name::Symbol)
    if name in fieldnames(Solution)
        return getproperty(getfield(estimator, :sol), name)
    else
        return getfield(estimator, name)
    end
end

# SOME INITIALIZATIONS

init_etacont(::OneD, ::Scalar, u::Vector{Float64}) = zero(u)