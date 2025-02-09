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
    @unpack dx, Nx = integrator.params.mesh
    @unpack u, uprev, dt, fnum, equation = integrator
    numflux!(integrator)
    # @views fluxforward = fnum[2:end,:]
    # @views fluxbackward = fnum[1:end-1,:]
    # @. u = uprev - dt / dx * (fluxforward - fluxbackward)

    for i in 1:Nx
        for j in 1:equation.p
            u[i,j] = uprev[i,j] - dt / dx * (fnum[i+1,j] - fnum[i,j])
        end
    end
end

function loopheader!(integrator::Integrator)
    dt_CFL!(integrator)
end

function loopfooter!(integrator::Integrator)
    integrator.t += integrator.dt
    integrator.niter += 1
    integrator.uprev .= integrator.u
    integrator.fcont .= integrator.equation.funcs.flux.(integrator.u)
    integrator.Dfcont .= integrator.equation.funcs.Dflux.(integrator.u)
    update_log!(integrator)
end

function update_log!(integrator::Integrator)
    @unpack log = integrator
    @unpack ulog, tlog, dtlog = log.config
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    ulog ? push!(log.ulog, copy(integrator.u)) : nothing
    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    tlog ? push!(log.tlog, integrator.t) : nothing 
    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dtlog ? push!(log.dtlog, integrator.dt) : nothing
end

