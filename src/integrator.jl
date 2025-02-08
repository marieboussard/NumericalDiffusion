struct IntegratorOptions
    maxiter::Int
end

mutable struct IntegratorCache <: Cache 

    cfl_loc::Float64
    sL::Int
    sR::Int
    stencil::Vector{Int}

    #IntegratorCache(sL, sR) = new(zero(Float64), sL, sR, zeros(Int, sL + sR))
    function IntegratorCache(sL, sR)
        new(zero(Float64), sL, sR, zeros(Int, sL + sR))
    end
end

mutable struct Integrator{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, scacheType <: Cache, tcacheTpe <: Cache, icacheType <: IntegratorCache}

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

    log::LogBook

    cfl::Float64
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config)
        
        # INIT SOLUTION AND FLUX
        uinit = equation.initcond.(params.mesh.x)
        #fnum = zeros(Float64, (params.mesh.Nx+1, equation.p))
        fnum = zeros(Float64, params.mesh.Nx+1)
        fcont = equation.flux.(uinit)
        uprev = copy(uinit)
        u = zero(uprev)
        
        # INIT SL AND SR
        sL, sR = compute_sL(time_scheme, space_scheme), compute_sR(time_scheme, space_scheme)

        # INTEGRATOR OPTIONS
        opts    = IntegratorOptions(maxiter)

        # INIT CACHE
        space_cache         = init_cache(space_scheme)
        time_cache          = init_cache(time_scheme)
        integrator_cache    = IntegratorCache(sL, sR)

        # INIT LOGBOOK
        logbook = LogBook(log_config)

        new{typeof(equation), typeof(params), typeof(time_scheme), typeof(space_scheme), typeof(u), typeof(space_cache), typeof(time_cache), typeof(integrator_cache)}(equation, params, time_scheme, space_scheme, u, uprev, uinit, fnum, fcont, 0, 0.0, params.t0, opts, space_cache, time_cache, integrator_cache, logbook, zero(Float64))
    end


end


