mutable struct Estimator

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
    method::QuantifMethod
    qparams::qparamsType

    # CACHE
    method_cache::mcacheType

    # RESULTS
    D::diffType

    function Estimator(sol::Solution, method::QuantifMethod, qparams::QuantifParams; kwargs...)
        entfun = entropy(typeof(sol.equation.eqfun))
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