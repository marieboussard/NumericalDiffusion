struct DiffEstimate{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, entfunType <: AbstractEntropyFun, etaType <: AbstractArray, dataType <: AbstractArray, methodType<:QuantifMethod, diffType<:AbstractArray, mcacheType<:MethodCache}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType
    space_scheme::sschemeType
    
    # ENTROPY
    entfun::entfunType
    etacont_init::etaType
    etacont::etaType

    # DATA
    uinit::dataType
    u::dataType
    dt::Float64                  # Final timestep
    t::Float64                   # Time reached

    # QUANTIFICATION METHOD
    method::methodType

    # DIFFUSION
    D::diffType

    # OTHER OUTPUTS
    output::mcacheType

    name::String
    function DiffEstimate(estimator::Estimator, name::String)
        new{typeof(estimator.equation), typeof(estimator.params), typeof(estimator.time_scheme), typeof(estimator.space_scheme), typeof(estimator.entfun), typeof(estimator.etacont), typeof(estimator.u), typeof(estimator.method), typeof(estimator.D), typeof(estimator.method_cache)}(estimator. equation, estimator.params, estimator.time_scheme, estimator.space_scheme, estimator.entfun, estimator.etacont_init, estimator.etacont, estimator.uinit, estimator.u, estimator.dt, estimator.t, estimator.method, estimator.D, estimator.method_cache, name)
    end
end

function Base.getproperty(estimate::DiffEstimate, name::Symbol)
    if name in fieldnames(typeof(getfield(estimate, :output)))
        return getproperty(getfield(estimate, :output), name)
    else
        return getfield(estimate, name)
    end
end