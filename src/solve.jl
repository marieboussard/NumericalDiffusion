function solve(equation, params, time_scheme, space_scheme; maxiter = 100, log_config = DefaultLogConfig, kwargs...)

    integrator = Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config; kwargs...)

    while integrator.t < integrator.params.tf && integrator.t < maxiter
        loopheader!(integrator)
        performstep!(integrator)
        loopfooter!(integrator)
    end
    Solution(integrator)
end

function performstep!(integrator::Integrator)
    @unpack dx = integrator.params.mesh
    @unpack u, uprev, dt, flux = integrator
    numflux!(integrator)
    @. u = uprev - dt / dx * (flux[2:end,:] .- flux[1:end-1,:])
end

function loopheader!(integrator::Integrator)
    @unpack problem = integrator
    @unpack params, scheme, equation = problem
    integrator.dt = min(compute_dt_with_CFL(integrator), params.Tf - integrator.t)
end

function loopfooter!(integrator::Integrator)
    integrator.t += dt
    integrator.niter += 1
    integrator.uprev .= integrator.u
    update_log!(integrator)
end

function update_log!(integrator::Integrator)
    @unpack log_book = integrator
    @unpack ulog, tlog, dtlog = log_book.config
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    ulog ? push!(log_book.u_log, copy(integrator.u)) : nothing
    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    tlog ? push!(log_book.t_log, integrator.t) : nothing 
    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dtlog ? push!(log_book.dt_log, integrator.dt) : nothing
end

