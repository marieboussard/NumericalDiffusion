struct IntegratorOptions
    maxiter::Int
end

struct IntegratorCache <: Cache end

mutable struct Integrator

    # PROBLEM COMPONENTS
    equation::Equation
    params::Parameters
    time_scheme::TimeScheme 
    space_scheme::SpaceScheme

    u
    uprev
    uinit
    flux

    niter::Int
    dt::Float64
    t::Float64

    opts::IntegratorOptions
    
    space_cache::Cache
    time_cache::Cache
    integrator_cache::Cache

    log::LogBook

    cfl::Float64
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config)
        uinit = equation.initcond.(params.mesh.x)
        uprev = copy(uinit)
        u = zero(uprev)

        flux = zeros(params.mesh.Nx+1, equation.p)

        new(equation, params, time_scheme, space_scheme, u, uprev, uinit, flux, 0, 0.0, params.t0, IntegratorOptions(maxiter), init_cache(space_scheme), init_cache(time_scheme), init_cache(), LogBook(log_config), zero(Float64))
    end

end


