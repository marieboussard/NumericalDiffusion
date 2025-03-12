mutable struct EstimatorCache
    sL::Int
    sR::Int
    indices
    utilde
    uhat
    function EstimatorCache(time_scheme::TimeScheme, space_scheme::SpaceScheme)
        sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)
        utilde = init_utilde()
        uhat = zero(utilde)
        new(sL, sR, utilde, uhat)
    end
end


mutable struct Estimator{entfunType<:AbstractEntropyFun, methodType<:QuantifMethod, mdtype<:ModifiedDataType, qparamsType <: QuantifParams, mcacheType<:MethodCache, gType <: AbstractArray}

    # # PROBLEM COMPONENTS
    # equation::equationType
    # params::parametersType
    # time_scheme::tschemeType 
    # space_scheme::sschemeType

    # # DATA
    # u0::dataType
    # u::dataType 
    # fnum::fnumType

    # # TIME
    # dt::Float64

    # PROBLEM COMPONENTS
    sol::Solution
    entfun::entfunType

    # QUANTIFICATION METHOD
    method::methodType
    md::mdtype
    qparams::qparamsType

    # CACHE
    cache::EstimatorCache
    method_cache::mcacheType

    # BOUNDS FOR NUMERICAL ENTROPY FLUX
    Gcont::gType
    m::gType
    M::gType

    # RESULTS
    D::diffType

    function Estimator(sol::Solution, method::QuantifMethod, qparams::QuantifParams; kwargs...)
        entfun = entropy(typeof(sol.equation.eqfun))

        # CACHE
        cache = EstimatorCache(sol.time_scheme, sol.space_scheme)
        method_cache = init_cache(method)

        new{}()
    end
end

# FASTER ACCESS TO SOME FIELDS OF ESTIMATOR

function Base.getproperty(estimator::Estimator, name::Symbol)
    if name in fieldnames(Solution)
        return getproperty(getfield(estimator, :sol), name)
    else
        return getfield(estimator, name)  # Accès normal aux autres champs de Outer
    end
end

# CACHE INITIALIZATION

init_cache(::Posteriori, args...) = PosterioriCache(args...)