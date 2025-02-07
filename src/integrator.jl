struct IntegratorOptions
    maxiter::Int
end

mutable struct IntegratorCache <: Cache 
    cfl_loc::Float64

    sL::Int
    sR::Int
    stencil

    IntegratorCache(sL, sR) = new(zero(Float64), sL, sR, zeros(Int, sL + sR))
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
    dt::Float64
    t::Float64

    opts::IntegratorOptions
    
    space_cache::Cache
    time_cache::Cache
    cache::Cache

    log::LogBook

    cfl::Float64
    

    function Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config)
        uinit = equation.initcond.(params.mesh.x)
        uprev = copy(uinit)
        u = zero(uprev)

        flux = zeros(Float64, (params.mesh.Nx+1, equation.p))
        sL, sR = compute_sL(time_scheme, space_scheme), compute_sR(time_scheme, space_scheme)

        new(equation, params, time_scheme, space_scheme, u, uprev, uinit, flux, 0, 0.0, params.t0, IntegratorOptions(maxiter), init_cache(space_scheme), init_cache(time_scheme), init_cache(sL, sR), LogBook(log_config), zero(Float64))
    end

end


