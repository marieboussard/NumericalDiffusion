struct DiffEstimate{methodType<:QuantifMethod, gType<:AbstractArray, diffType<:AbstractArray, mcacheType<:MethodCache}

    # # PROBLEM COMPONENTS
    # equation::Equation
    # params::Parameters
    # time_scheme::TimeScheme
    # space_scheme::SpaceScheme

    # # DATA
    # uinit::utype
    # u::utype

    # QUANTIFICATION METHOD
    method::methodType

    # BOUNDS FOR NUMERICAL ENTROPY FLUX
    m::gType
    M::gType

    # DIFFUSION
    D::diffType

    # OTHER OUTPUTS
    output::mcacheType

    name::String
    function DiffEstimate(estimator::Estimator, name::String)
        new{typeof(estimator.method), typeof(estimator.m), typeof(estimator.D), typeof(estimator.method_cache)}(estimator.method, estimator.m, estimator.M, estimator.D, estimator.method_cache, name)
    end
end

function Base.getproperty(estimate::DiffEstimate, name::Symbol)
    if name in fieldnames(typeof(getfield(estimate, :output)))
        return getproperty(getfield(estimate, :output), name)
    else
        return getfield(estimate, name)
    end
end