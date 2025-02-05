function solve(equation, params, timeScheme, spaceScheme; NiterMax=100, kwargs...)
    integrator = Integrator()
    while integrator.sol.t < problem.params.tf && integrator.sol < NiterMax
        loopheader!(integrator)
        performstep!(integrator)
        loopfooter!(integrator)
    end
    integrator.sol
end

solve(problem) = solve(problem.equation, problem.params, problem.timeScheme, problem.spaceScheme)

function initialize_FV(problem::FVProblem)
    @unpack domain, equation, scheme, saveLog, u_init = problem
    if saveLog
        FVSolution(problem, copy(u_init), 0, domain.t0, 0.0, FVLog([copy(u_init)], typeof(t0)[], [t0]))
    else
        FVSolution(problem, copy(u_init), 0, domain.t0, 0.0, nothing)
    end
end

function performstep!(fv_sol::FVSolution)
    @unpack problem, u_approx, Nt, t = fv_sol
    @unpack domain, equation, scheme = problem
    @unpack dx, Tf = domain
    # Find the next time step with a CFL condition
    
    numericalFluxMat = vecNumFlux(equation.source, scheme, equation, collect(u_approx); dt=dt, domain=domain)
    u_approx .= u_approx .- dt / domain.dx * (numericalFluxMat[2:end,:] .- numericalFluxMat[1:end-1,:])
    
end

function loopheader!(integrator::Integrator)
    @unpack u = integrator
    dt = min(compute_dt_with_CFL(scheme, equation, u, dx), Tf - t)
end

function loopfooter!(integrator::Integrator)
    fv_sol.u_approx = u_approx
    fv_sol.dt = dt
    fv_sol.t += dt
    fv_sol.Nt += 1

    save_step!(integrator)
end

function save_step!(integrator::Integrator)
    @unpack log = fv_sol
    if isnothing(log)
        @error "Trying to update log, but it has not been initialized."
    end
    @unpack u_log, dt_log, t_log = log
    @unpack u_approx, dt, t = fv_sol
    push!(u_log, copy(u_approx))
    push!(dt_log, dt)
    push!(t_log, t)
end

