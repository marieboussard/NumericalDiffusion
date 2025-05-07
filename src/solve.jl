function solve(equation::Equation, params::Parameters, time_scheme::TimeScheme, space_scheme::SpaceScheme; maxiter::Int = 10000, log_config::LogConfig = DefaultLogConfig, name::String="", kwargs...)

    integrator = Integrator(equation, params, time_scheme, space_scheme, maxiter, log_config; kwargs...)
    initialize_integrator!(integrator)

    while integrator.t < integrator.params.tf && integrator.niter < maxiter
        loopheader!(integrator)
        performstep!(integrator.equation.dim, integrator)
        loopfooter!(integrator)
    end
    Solution(integrator, name)
end

function performstep!(::OneD, integrator::Integrator)
    @unpack dx, Nx = integrator.params.mesh
    @unpack u, uprev, dt, fnum, equation, time_scheme = integrator
    
    # numflux!(integrator)
    numflux!(time_scheme, integrator)
    for i in 1:Nx
        for j in 1:equation.p
            u[i,j] = uprev[i,j] - dt / dx * (fnum[i,j] - fnum[mod1(i-1, Nx),j])
        end
    end

    if has_source(equation.source)
        discretize_sourceterm!(equation.dim, equation.source.source_discretize, integrator)
        for i in 1:Nx
            for j in 1:equation.p
                u[i,j] += dt * integrator.cache.sourceterm[i,j]
            end
        end
    end
end

function performstep!(::TwoD, integrator::Integrator)
    @unpack dx, dy, Nx, Ny = integrator.params.mesh
    @unpack u, uprev, dt, fnum, equation = integrator
    numflux2D!(integrator)

    for j in 1:Nx
        for k in 1:Ny
            for r in 1:equation.p
            u[j,k,r] = uprev[j,k,r] - dt/dx *(fnum[j,k,r,1]-fnum[mod1(j-1,Nx),k,r,1]) - dt/dy * (fnum[j,k,r,2]-fnum[j,mod1(k-1,Nx),r,2])
            end
        end
    end

    if has_source(equation.source)
        discretize_sourceterm!(equation.dim, equation.source.source_discretize, integrator)
        for j in 1:Nx
            for k in 1:Ny
                for r in 1:equation.p
                    u[j,k,r] += dt*integrator.cache.sourceterm[j,k,r]
                end
            end
        end
    end
end


function loopheader!(integrator::Integrator)
    dt_CFL!(integrator.equation.dim, integrator)
end

function loopfooter!(integrator::Integrator)
    @unpack equation, uprev, u = integrator
    @unpack source = equation
    integrator.t += integrator.dt
    integrator.niter += 1
    uprev .= u
    update_log!(integrator)
    flux!(equation.dim, equation.funcs, integrator)
    update_cflcache!(equation.dim, equation.eqtype, equation.funcs, integrator)
    
    nothing
end

# function update_flux!(::OneD, integrator::Integrator)
#     @unpack equation, u = integrator
#     @views vu = u
#     # fcont .= flux(equation.funcs, vu)
#     flux!(equation.funcs, vu, integrator.fcont)
# end

function update_cflcache!(::OneD, ::Scalar, eqfun::AbstractEquationFun, integrator)
    Dflux!(eqfun, integrator.u, integrator.cache.cfl_cache.absDfcont)
    abs!(integrator.cache.cfl_cache.absDfcont, integrator.cache.cfl_cache.absDfcont)
end

function update_cflcache!(::OneD, ::Scalar, eqfun::AbstractEquationFun, u::AbstractVector, cfl_cache::CFLCache)
    @unpack absDfcont = cfl_cache
    Dflux!(eqfun, u, absDfcont)
    abs!(absDfcont, absDfcont)
end

function update_cflcache!(::TwoD, ::Scalar, eqfun::AbstractEquationFun, integrator)
    @unpack Dfcont = integrator.cache.cfl_cache
    Dflux!(eqfun, integrator.u, selectdim(Dfcont, ndims(Dfcont),1), selectdim(Dfcont, ndims(Dfcont),2))
end

# update_cflcache!(eqtype::EquationType, eqfun::AbstractEquationFun, integrator::Integrator) = nothing

function update_log!(integrator::Integrator)
    @unpack log = integrator
    @unpack ulog, tlog, dtlog, fnumlog, fcontlog = log.config
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    ulog ? push!(log.ulog, copy(integrator.u)) : nothing
    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    tlog ? push!(log.tlog, integrator.t) : nothing 
    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dtlog ? push!(log.dtlog, integrator.dt) : nothing
    # STORE INTERMEDIATE NUMERICAL FLUX
    fnumlog ? push!(log.fnumlog, integrator.fnum) : nothing
    # STORE INTERMEDIATE POINTWISE FLUX FUNCTION
    fcontlog ? push!(log.fcontlog, integrator.fcont) : nothing
    nothing
end

