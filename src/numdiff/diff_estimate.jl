struct DiffEstimate{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, methodType<:QuantifMethod, gType<:AbstractArray, diffType<:AbstractArray, mcacheType<:MethodCache}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType
    space_scheme::sschemeType

    # DATA
    uinit::dataType
    u::dataType
    dt::Float64                  # Final timestep
    t::Float64                   # Time reached

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
        new{typeof(estimator.equation), typeof(estimator.params), typeof(estimator.time_scheme), typeof(estimator.space_scheme), typeof(estimator.u), typeof(estimator.method), typeof(estimator.m), typeof(estimator.D), typeof(estimator.method_cache)}(estimator. equation, estimator.params, estimator.time_scheme, estimator.space_scheme, estimator.uinit, estimator.u, estimator.dt, estimator.t, estimator.method, estimator.m, estimator.M, estimator.D, estimator.method_cache, name)
    end
end

function Base.getproperty(estimate::DiffEstimate, name::Symbol)
    if name in fieldnames(typeof(getfield(estimate, :output)))
        return getproperty(getfield(estimate, :output), name)
    else
        return getfield(estimate, name)
    end
end