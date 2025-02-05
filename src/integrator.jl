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

end


