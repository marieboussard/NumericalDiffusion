struct IntegratorOptions

    NiterMax::Int

end

mutable struct Integrator

    prob::Problem
    u
    u_prev
    flux

    niter::Int
    dt
    t

    opts::IntegratorOptions
    
    space_cache::Cache
    time_cache::Cache

    log_config::LogConfig
    log_book::LogBook
    

    function Integrator(problem, opts, log_config)
        @unpack equation, params = problem
        u_prev, u = copy(equation.u_init), copy(equation.u_init)

        flux = zeros(params.mesh.Nx, equation.p)

        new(problem, u, u_prev, flux, 0, 0.0, params.t0, opts, initialize_cache(problem.spaceScheme), initialize_cache(problem.timeScheme), log_config, LogBook(log_config))
    end

end


