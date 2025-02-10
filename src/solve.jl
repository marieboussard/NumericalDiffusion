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
    @unpack u, uprev, dt, fnum, equation, sourceterm = integrator
    numflux!(integrator)
    # @views fluxforward = fnum[2:end,:]
    # @views fluxbackward = fnum[1:end-1,:]
    # @. u = uprev - dt / dx * (fluxforward - fluxbackward)

    for i in 1:Nx
        for j in 1:equation.p
            u[i,j] = uprev[i,j] - dt / dx * (fnum[i+1,j] - fnum[i,j])
        end
    end

    if has_source(equation.source)
        for i in 1:Nx
            for j in 1:equation.p
                u[i,j] += dt * integrator.sourceterm[i,j]
            end
        end
    end
end

function loopheader!(integrator::Integrator)
    dt_CFL!(integrator)
end

function loopfooter!(integrator::Integrator)
    @unpack cache, equation, uprev, u, fcont = integrator
    @unpack source = equation

    integrator.t += integrator.dt
    integrator.niter += 1
    uprev .= u
    fcont .= flux(equation.funcs, u)
    # STORING FLUX DERIVATIVE IN SCALAR CASE
    if equation.p ==1
        cache.Dfcont .= Dflux(equation.funcs, u)
    end
    # UPDATING SOURCE TERM
    if has_source(source)
        discretize_sourceterm!(source.source_discretize, integrator)
    end
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

