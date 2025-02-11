struct IntegratorOptions
    maxiter::Int
end

mutable struct IntegratorCache{typeDfcont, sourcetermType} <: Cache 

    cfl_loc::Float64
    sL::Int
    sR::Int
    stencil::Vector{Int}
    Dfcont::typeDfcont
    sourceterm::sourcetermType

    function IntegratorCache(sL, sR, equation, uinit, x)
        Dfcont = init_Dfcont(equation.eqtype, equation, uinit)
        sourceterm = init_sourceterm(equation.source, uinit, x)
        new{typeof(Dfcont), typeof(sourceterm)}(zero(Float64), sL, sR, zeros(Int, sL + sR), Dfcont, sourceterm)
    end
end

# INIT CACHE CONTENT
init_Dfcont(::Scalar, equation, uinit) = Dflux(equation.funcs, uinit)
init_Dfcont(::System, equation, uinit) = nothing
init_sourceterm(::NoSource, args...) = nothing


mutable struct Integrator{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, scacheType <: Cache, tcacheTpe <: Cache, icacheType <: IntegratorCache, srcacheType<:sourceCacheType}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType 
    space_scheme::sschemeType

    u::dataType
    uprev::dataType
    uinit::dataType
    fnum::dataType
    fcont::dataType

    niter::Int
    dt::Float64
    t::Float64

    opts::IntegratorOptions
    
    space_cache::scacheType
    time_cache::tcacheType
    cache::icacheType
    source_cache::srcacheType

    log::LogBook

    cfl::Float64
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config::LogConfig)
        
        # INIT SOLUTION AND FLUX
        uinit = equation.initcond(params.mesh.x)
        if equation.p == 1
            fnum = zeros(Float64, params.mesh.Nx+1)
        else
            fnum = zeros(Float64, (params.mesh.Nx+1, equation.p))
        end
        fcont = flux(equation.funcs, uinit)
        uprev = copy(uinit)
        u = zero(uprev)
        
        # INIT SL AND SR
        sL, sR = compute_sL(time_scheme, space_scheme), compute_sR(time_scheme, space_scheme)

        # INTEGRATOR OPTIONS
        opts    = IntegratorOptions(maxiter)

        # INIT CACHE
        space_cache         = init_cache(space_scheme)
        time_cache          = init_cache(time_scheme)
        integrator_cache    = IntegratorCache(sL, sR, equation, uinit, params.mesh.x)
        source_cache        = init_cache(equation.source, equation.source.source_discretize, params.mesh.x)

        # INIT LOGBOOK
        logbook = LogBook(log_config)

        new{typeof(equation), typeof(params), typeof(time_scheme), typeof(space_scheme), typeof(u), typeof(space_cache), typeof(time_cache), typeof(integrator_cache), typeof(source_cache)}(equation, params, time_scheme, space_scheme, u, uprev, uinit, fnum, fcont, 0, 0.0, params.t0, opts, space_cache, time_cache, integrator_cache, source_cache, logbook, zero(Float64))
    end


end


# # INIT INTEGRATOR CONTENT
# init_source_discretize(::NoSource) = nothing
# init_source_discretize()


# function Integrator(equation::equationType, params::paramsType, time_scheme::tschemeType, space_scheme::sschemeType, maxiter::Int, log_config::LogConfig) where equationType<:Equation where paramsType<:Parameters where tschemeType<:TimeScheme where sschemeType<:SpaceScheme
#         Integrator(equation, params, time_scheme, space_scheme, maxiter, nothing, log_config)
# end