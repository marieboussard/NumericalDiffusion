struct IntegratorOptions

    maxiter::Int

end

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
    dt
    t

    opts::IntegratorOptions
    
    space_cache::Cache
    time_cache::Cache
    integrator_cache::Cache

    log::LogBook
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config)
        uinit = equation.initcond.(params.mesh.x)
        u_prev = copy(uinit)
        u = zero(u_prev)

        flux = zeros(params.mesh.Nx, equation.p)

        new(equation, params, time_scheme, space_scheme, u, uprev, uinit, flux, 0, 0.0, params.t0, IntegratorOptions(maxiter), initialize_cache(space_scheme), init_cache(time_scheme), init_cache(integrator), LogBook(log_config))
    end

end


