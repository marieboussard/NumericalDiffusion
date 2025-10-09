struct IntegratorOptions
    maxiter::Int
end

mutable struct IntegratorCache{cflCacheType<:CFLCache, sourcetermType} <: Cache 

    sL::Int
    sR::Int
    #stencil::Vector{Int}
    cfl_cache::cflCacheType
    sourceterm::sourcetermType

    function IntegratorCache(sL::Int, sR::Int, equation::Equation, uinit::AbstractArray, mesh::Mesh, source_cache)
        # Dfcont = init_Dfcont(equation.eqtype, equation, uinit)
        cfl_cache = init_cfl_cache(equation.dim, equation.eqtype, equation.funcs, equation, uinit)
        #sourceterm = init_sourceterm(equation.source, uinit, mesh, source_cache)
        sourceterm = init_sourceterm(equation.source, uinit)
        # new{typeof(cfl_cache), typeof(sourceterm)}(sL, sR, zeros(Int, sL + sR), cfl_cache, sourceterm)
        new{typeof(cfl_cache), typeof(sourceterm)}(sL, sR, cfl_cache, sourceterm)
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
init_sourceterm(::AbstractSource, uinit, args...) = zero(uinit)

"Contain all information necessary for finite volume solving: equation, parameters, scheme, initial data, etc."
mutable struct Integrator{equationType <: Equation, parametersType <: Parameters, tschemeType <: TimeScheme, sschemeType <: SpaceScheme, dataType <: AbstractArray, fnumType<:AbstractArray, tcacheType <: TimeCache, srcacheType, icacheType <: IntegratorCache}

    # PROBLEM COMPONENTS
    equation::equationType
    params::parametersType
    time_scheme::tschemeType 
    space_scheme::sschemeType

    # DATA
    u::dataType
    uprev::dataType
    uinit::dataType
    fnum::fnumType
    fcont::fnumType

    # TIME PARAMETERS
    niter::Int
    dt::Float64
    t::Float64

    # OPTIONS
    opts::IntegratorOptions
    
    # CACHE
    # space_cache::scacheType
    time_cache::tcacheType
    source_cache::srcacheType
    cache::icacheType

    # LOGBOOK
    log::LogBook

    function Integrator(equation::Equation, params::Parameters, time_scheme::TimeScheme, space_scheme::SpaceScheme, maxiter::Int, log_config::LogConfig)
        
        # INIT SOLUTION AND FLUX
        uinit = initialize_u(equation.dim, equation.eqtype, equation.source, equation, params)
        fnum = init_flux(equation.dim, equation.eqtype, equation, params.mesh)
        fcont = zero(fnum)
        uprev = copy(uinit)
        #u = zero(uprev)
        u = copy(uprev)
        
        # INIT SL AND SR
        sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)

        # INTEGRATOR OPTIONS
        opts    = IntegratorOptions(maxiter)

        # INIT CACHE
        #space_cache         = init_cache(space_scheme, params.mesh.Nx, equation.dim)
        time_cache          = init_cache(time_scheme, space_scheme, uinit, equation.dim, params, 0.0)
        source_cache        = init_cache(equation.source, params.mesh)
        integrator_cache    = IntegratorCache(sL, sR, equation, uinit, params.mesh, source_cache)

        # INIT LOGBOOK
        logbook = LogBook(log_config, u, fnum)

        new{typeof(equation), typeof(params), typeof(time_scheme), typeof(space_scheme), typeof(u), typeof(fnum), typeof(time_cache), typeof(source_cache), typeof(integrator_cache)}(equation, params, time_scheme, space_scheme, u, uprev, uinit, fnum, fcont, 0, 0.0, params.t0, opts, time_cache, source_cache, integrator_cache, logbook)
    end


end


# # INIT INTEGRATOR CONTENT

initialize_u(::OneD, ::Scalar, ::NoSource, equation::AbstractEquation, params::Parameters, args...) = equation.initcond(params.mesh.x)

function initialize_u(::OneD, ::System, ::NoSource, equation::AbstractEquation, params::Parameters, args...)
    uinit = zeros(eltype(params.mesh.x), (params.mesh.Nx, equation.p))
    for j in 1:params.mesh.Nx
        uinit[j,:] .= equation.initcond(params.mesh.x[j])
    end
    uinit
end

function initialize_u(::TwoD, ::Scalar, ::NoSource, equation::AbstractEquation, params::Parameters)
    @unpack Nx, Ny, x, y = params.mesh
    uinit = zeros(eltype(x), (Nx, Ny))
    for j in eachindex(x)
        for k in eachindex(y)
            uinit[j,k] = equation.initcond(x[j], y[k])
        end
    end
    uinit
end

function initialize_u(::TwoD, ::System, ::NoSource, equation::AbstractEquation, params::Parameters)
    @unpack Nx, Ny, x, y = params.mesh
    uinit = zeros(eltype(x), (Nx, Ny, equation.p))
    for j in eachindex(x)
        for k in eachindex(y)
            uinit[j,k,:] .= equation.initcond(x[j], y[k])
        end
    end
    uinit
end

# init_flux(::OneD, ::Scalar, equation::Equation, mesh::OneDMesh) = zeros(Float64, mesh.Nx+1)
# init_flux(::OneD, ::System, equation::Equation, mesh::OneDMesh) = zeros(Float64, (mesh.Nx+1, equation.p))
init_flux(::OneD, ::Scalar, equation::Equation, mesh::OneDMesh) = zeros(Float64, mesh.Nx)
init_flux(::OneD, ::System, equation::Equation, mesh::OneDMesh) = zeros(Float64, (mesh.Nx, equation.p))
init_flux(::TwoD, ::EquationType, equation::Equation, mesh::Mesh) = zeros(Float64, mesh.Nx, mesh.Ny, equation.p, 2)

# FILL INTEGRATOR FIELDS WITH INITIAL VALUES

function initialize_integrator!(integrator::Integrator)
    @unpack equation = integrator
    flux!(equation.dim, equation.funcs, integrator)
    fillcache!(integrator.space_scheme, integrator)
    nothing
end