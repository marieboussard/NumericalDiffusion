mutable struct EstimatorCache{utype <: AbstractArray, cflCacheType<:CFLCache, sourcetermType, mdcacheType<:ModifiedDataCache} <: Cache
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
    sourceterm_tilde::sourcetermType
    mdcache::mdcacheType
    
    function EstimatorCache(equation::Equation, time_scheme::TimeScheme, space_scheme::SpaceScheme, u::AbstractArray, method::QuantifMethod)
        sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)
        indices = init_indices(method.mdtype, method.boundstype, sL, sR)
        utilde = init_utilde(method.mdtype, method.boundstype, equation.dim, equation.eqtype, u, sL, sR)
        uhat = init_uhat(method.mdtype, method.boundstype, equation.dim, equation.eqtype, u, sL, sR)
        fcont_tilde = zero(utilde)
        ftilde = init_ftilde(method.mdtype, method.boundstype, equation.dim, equation.eqtype, u, sL, sR)
        cfl_cache = init_cfl_cache(equation.dim, equation.eqtype, equation.funcs, equation, utilde)
        eta_tilde = zero(uhat)
        eta_hat = zero(uhat)
        sourceterm_tilde = init_sourceterm(equation.source, utilde)
        mdcache = init_cache(method.mdtype, equation, u)

        new{typeof(utilde), typeof(cfl_cache), typeof(sourceterm_tilde), typeof(mdcache)}(sL, sR, indices, utilde, uhat, fcont_tilde, ftilde, cfl_cache, eta_tilde, eta_hat, sourceterm_tilde, mdcache)
    end
end




mutable struct Estimator{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, methodType<:QuantifMethod, ecacheType<:EstimatorCache, mcacheType<:MethodCache, scacheType<:SpaceCache, tcacheType<:TimeCache, srcacheType, entfunType<:AbstractEntropyFun, etaType <: AbstractArray, diffType<:AbstractArray}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType
    space_scheme::sschemeType

    # sol::Solution

    # DATA
    uinit::dataType
    u::dataType
    dt::Float64                  # Final timestep
    t::Float64                   # Time reached
    
    # QUANTIFICATION METHOD
    method::methodType

    # CACHE
    cache::ecacheType
    method_cache::mcacheType
    space_cache::scacheType
    time_cache::tcacheType
    source_cache::srcacheType
    
    # ENTROPY
    entfun::entfunType
    etacont_init::etaType
    etacont::etaType

    # RESULTS
    D::diffType

    function Estimator(sol::Solution, method::QuantifMethod; kwargs...)

        # INIT PROBLEM COMPONENTS
        equation = sol.equation
        params = sol.params
        time_scheme = sol.time_scheme
        space_scheme = sol.space_scheme

        # INIT DATA
        if sol.niter == 1
            uinit = sol.uinit
            dt = sol.dt
        elseif sol.log.config.ulog && sol.log.config.dtlog
            uinit = sol.log.ulog[end-1]
            dt = sol.log.dtlog[end]
        else
            throw("Diffusion quantification can only be done for a single timestep : please give a solution with a single iteration or with a saving of intermediate values")
        end
        u = sol.u
        t = sol.t

        # INIT CACHE
        cache = EstimatorCache(sol.equation, sol.time_scheme, sol.space_scheme, sol.u, method)
        method_cache = init_cache(method, sol.equation, sol.u)
        space_cache = init_cache(sol.space_scheme)
        time_cache = init_cache(sol.time_scheme)
        source_cache = init_cache(sol.equation.source, sol.params.mesh)
        
        # INIT ENTROPY
        entfun = entropy(typeof(sol.equation.funcs))
        etacont_init = init_etacont(sol.equation.dim, sol.equation.eqtype, sol.u)
        etacont = zero(etacont_init)

        # INIT DIFFUSION
        D = zero(etacont_init)

        new{typeof(equation), typeof(params), typeof(time_scheme), typeof(space_scheme), typeof(u), typeof(method), typeof(cache), typeof(method_cache), typeof(space_cache), typeof(time_cache), typeof(source_cache), typeof(entfun), typeof(etacont), typeof(D)}(equation, params, time_scheme, space_scheme, uinit, u, dt, t, method, cache, method_cache, space_cache, time_cache, source_cache, entfun, etacont_init, etacont, D)
    end
end

# # FASTER ACCESS TO SOME FIELDS OF ESTIMATOR

# function Base.getproperty(estimator::Estimator, name::Symbol)
#     if name in fieldnames(Solution)
#         return getproperty(getfield(estimator, :sol), name)
#     else
#         return getfield(estimator, name)
#     end
# end

# SOME INITIALIZATIONS

init_etacont(::OneD, ::Scalar, u::Vector{Float64}) = zero(u)
init_etacont(::OneD, ::System, u::Matrix{Float64}) = zero(selectdim(u,2,1))

# FILL ESTIMATOR FIELDS WITH INITIAL VALUES

function initialize_estimator!(estimator::Estimator)
    eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
    eta!(estimator.entfun, estimator.u, estimator.etacont)
    if has_source(estimator.equation.source)
        etasource!(estimator.entfun, estimator.uinit, estimator.source_cache, estimator.etacont_init)
        etasource!(estimator.entfun, estimator.u, estimator.source_cache, estimator.etacont)
    end
end