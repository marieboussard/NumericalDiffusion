function solve(equation, params, time_scheme, space_scheme; maxiter = 100, log_config = DefaultLogConfig, name="", kwargs...)

    integrator = Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config; kwargs...)

    while integrator.t < integrator.params.tf && integrator.t < maxiter
        loopheader!(integrator)
        performstep!(integrator)
        loopfooter!(integrator)
    end
    Solution(integrator, name)
end

function performstep!(integrator::Integrator)
    @unpack dx = integrator.params.mesh
    @unpack u, uprev, dt, fnum = integrator
    numflux!(integrator.time_scheme, integrator)
    @views fluxforward = fnum[2:end,:]
    @views fluxbackward = fnum[1:end-1,:]
    @. u = uprev - dt / dx * (fluxforward .- fluxbackward)
end

function loopheader!(integrator::Integrator)
    dt_CFL!(integrator)
end

function loopfooter!(integrator::Integrator)
    integrator.t += integrator.dt
    integrator.niter += 1
    integrator.uprev .= integrator.u
    integrator.fcont .= integrator.equation.flux.(integrator.u)
    update_log!(integrator)
end

function update_log!(integrator::Integrator)
    @unpack log = integrator
    @unpack ulog, tlog, dtlog = log.config
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    ulog ? push!(log.u_log, copy(integrator.u)) : nothing
    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    tlog ? push!(log.t_log, integrator.t) : nothing 
    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dtlog ? push!(log.dt_log, integrator.dt) : nothing
end

