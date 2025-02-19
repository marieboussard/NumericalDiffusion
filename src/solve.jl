function solve(equation, params, time_scheme, space_scheme; maxiter = 1000, log_config = DefaultLogConfig, name="", kwargs...)

    integrator = Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config; kwargs...)

    while integrator.t < integrator.params.tf && integrator.niter < maxiter
        loopheader!(integrator)
        performstep!(integrator.equation.dim, integrator)
        loopfooter!(integrator)
    end
    Solution(integrator, name)
end

function performstep!(::OneD, integrator::Integrator)
    @unpack dx, Nx = integrator.params.mesh
    @unpack u, uprev, dt, fnum, equation, cache = integrator
    @unpack sourceterm = cache
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
                u[i,j] += dt * sourceterm[i,j]
            end
        end
    end
end

function performstep!(::TwoD, integrator::Integrator)
    @unpack dx, Nx, Ny = integrator.params.mesh
    @unpack u, uprev, dt, fnum, equation, cache = integrator
    @unpack sourceterm = cache
    numflux!(integrator)

    for j in 1:Nx
        for k in 1:Ny
            for r in 1:equation.p
            u[j,k,r] = uprev[j,k,r] - dt/dx *(fnum[j+1,k,r]-fnum[j,k,r]) - dt/dy * (hnum[j,k+1,r]-hnum[j,k,r])
            end
        end
    end

    if has_source(equation.source)
        for j in 1:Nx
            for k in 1:Ny
                for r in 1:equation.p
                    u[j,k,r] += dt*sourceterm[j,k,r]
                end
            end
        end
    end
end


function loopheader!(integrator::Integrator)
    dt_CFL!(integrator)
end

function loopfooter!(integrator::Integrator)
    @unpack equation, uprev, u, fcont = integrator
    @unpack source = equation

    integrator.t += integrator.dt
    integrator.niter += 1
    uprev .= u
    fcont .= flux(equation.funcs, u)
    # STORING FLUX DERIVATIVE IN SCALAR CASE
    # if equation.p ==1
    #     cache.cfl_cache.Dfcont .= Dflux(equation.funcs, u)
    # end
    update_cflcache!(equation.eqtype, equation.funcs, integrator)
    # UPDATING SOURCE TERM
    if has_source(source)
        discretize_sourceterm!(source.source_discretize, integrator)
    end
    update_log!(integrator)
end

function update_cflcache!(::Scalar, eqfun::AbstractEquationFun, integrator)
    integrator.cache.cfl_cache.Dfcont .= Dflux(eqfun, integrator.u)
end

# update_cflcache!(eqtype::EquationType, eqfun::AbstractEquationFun, integrator::Integrator) = nothing

function update_log!(integrator::Integrator)
    @unpack log = integrator
    @unpack ulog, tlog, dtlog, fnumlog = log.config
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    ulog ? push!(log.ulog, copy(integrator.u)) : nothing
    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    tlog ? push!(log.tlog, integrator.t) : nothing 
    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dtlog ? push!(log.dtlog, integrator.dt) : nothing
    # STORE INTERMEDIATE NUMERICAL FLUX
    fnumlog ? push!(log.fnumlog, integrator.fnum) : nothing
end

