struct IntegratorOptions
    maxiter::Int
end

mutable struct IntegratorCache{cflCacheType<:CFLCacheType, sourcetermType} <: Cache 

    sL::Int
    sR::Int
    stencil::Vector{Int}
    cfl_cache::cflCacheType
    sourceterm::sourcetermType

    function IntegratorCache(sL, sR, equation, uinit, x)
        # Dfcont = init_Dfcont(equation.eqtype, equation, uinit)
        cfl_cache = init_cfl_cache(equation.dim, equation.eqtype, equation.funcs, equation, uinit)
        sourceterm = init_sourceterm(equation.source, uinit, x)
        new{typeof(cfl_cache), typeof(sourceterm)}(sL, sR, zeros(Int, sL + sR), cfl_cache, sourceterm)
    end
end

# INIT CACHE CONTENT
# function init_Dfcont(::Scalar, equation, uinit)
#     Dfcont= zero(uinit)
#     Dfcont .= Dflux(equation.funcs, uinit)
# end
# init_Dfcont(::System, equation, uinit) = nothing

init_cfl_cache(::OneD, ::Scalar, ::AbstractEquationFun, args...) = CFLCacheScalar(args...)
init_cfl_cache(::TwoD, ::Scalar, ::AbstractEquationFun, args...) = CFLCacheScalar2D(args...)
# init_cfl_cache(::System, equation, args...) = init_cfl_cache(equation, args...)
init_sourceterm(::NoSource, args...) = nothing
# init_sourceterm(source::AbstractSource, args...) = init_sourceterm(source, source.source_discretize, args...)


mutable struct Integrator{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, fnumType, fcontType, scacheType <: Cache, tcacheTpe <: Cache, icacheType <: IntegratorCache, srcacheType}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType 
    space_scheme::sschemeType

    u::dataType
    uprev::dataType
    uinit::dataType
    fnum::fnumType
    fcont::fcontType

    niter::Int
    dt::Float64
    t::Float64

    opts::IntegratorOptions
    
    space_cache::scacheType
    time_cache::tcacheType
    cache::icacheType
    source_cache::srcacheType

    log::LogBook

    # cfl::Float64
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config::LogConfig)
        
        # INIT SOLUTION AND FLUX
        uinit = initialize_u(equation.dim, equation.eqtype, equation.source, equation, params)
        fnum = init_fnum(equation.dim, equation, params.mesh)
        fcont = init_fcont(equation.dim, equation, uinit)
        uprev = copy(uinit)
        u = zero(uprev)
        
        # INIT SL AND SR
        sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)

        # INTEGRATOR OPTIONS
        opts    = IntegratorOptions(maxiter)

        # INIT CACHE
        space_cache         = init_cache(space_scheme)
        time_cache          = init_cache(time_scheme)
        integrator_cache    = IntegratorCache(sL, sR, equation, uinit, params.mesh.x)
        source_cache        = init_cache(equation.source, params.mesh.x)

        # INIT LOGBOOK
        logbook = LogBook(log_config)

        new{typeof(equation), typeof(params), typeof(time_scheme), typeof(space_scheme), typeof(u), typeof(fnum), typeof(fcont), typeof(space_cache), typeof(time_cache), typeof(integrator_cache), typeof(source_cache)}(equation, params, time_scheme, space_scheme, u, uprev, uinit, fnum, fcont, 0, 0.0, params.t0, opts, space_cache, time_cache, integrator_cache, source_cache, logbook)
    end


end


# # INIT INTEGRATOR CONTENT

initialize_u(::OneD, ::Scalar, ::NoSource, equation::AbstractEquation, params::Parameters, args...) = equation.initcond(params.mesh.x)
function initialize_u(::OneD, ::System,  ::NoSource, equation::AbstractEquation, params::Parameters, args...)
    uinit = zeros(eltype(x), (params.mesh.Nx, equation.p))
    for j in 1:Nx
        uinit[j,:] .= equation.initcond(params.mesh.x)
    end
end
init_fnum(::OneD, args...) = OneDFnum(args...).fnum
init_fnum(::TwoD, args...) = TwoDFnum(args...)
init_fcont(::OneD, args...) = OneDFcont(args...).fcont
init_fcont(::TwoD, args...) = TwoDFcont(args...)

# init_source_discretize(::NoSource) = nothing
# init_source_discretize()


# function Integrator(equation::equationType, params::paramsType, time_scheme::tschemeType, space_scheme::sschemeType, maxiter::Int, log_config::LogConfig) where equationType<:Equation where paramsType<:Parameters where tschemeType<:TimeScheme where sschemeType<:SpaceScheme
#         Integrator(equation, params, time_scheme, space_scheme, maxiter, nothing, log_config)
# end


# function initialize_integrator(integrator::Integrator)
#     update_flux!(equation.dim, integrator)
#     update_cflcache!(equation.dim, equation.eqtype, equation.funcs, integrator)
# end