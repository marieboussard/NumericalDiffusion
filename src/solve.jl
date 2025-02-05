function solve(problem; maxiter=100, name="", log_config=DefaultLogConfig, kwargs...)

    opts = IntegratorOptions(maxiter)
    integrator = Integrator(problem, opts, log_config)

    @unpack tf = problem.params
     
    while integrator.t < tf && integrator.t < maxiter
        loopheader!(integrator)
        performstep!(integrator)
        loopfooter!(integrator)
    end
    Solution(problem, status, integrator.niter, integrator.u, integrator.dt, integrator.t, integrator.log_book, name)
end

solve(equation, params, timeScheme, spaceScheme; kwargs...) = solve(Problem(equation, params, timeScheme, spaceScheme); kwargs...)

function performstep!(integrator::Integrator)
    @unpack dx = integrator.problem.params.mesh
    
    integrator.flux = vec_num_flux(integrator)
    integrator.u .= u .- integrator.dt / dx * (integrator.flux[2:end,:] .- integrator.flux[1:end-1,:])
    
end


function performstep_2D!(integrator::Integrator)

end



function loopheader!(integrator::Integrator)
    @unpack problem = integrator
    @unpack params = problem
    integrator.dt = min(compute_dt_with_CFL(problem.scheme, problem.equation, integrator.u_prev, params.mesh.dx), params.Tf - integrator.t)
end

function loopfooter!(integrator::Integrator)
    integrator.t += dt
    integrator.niter += 1
    integrator.u_prev .= integrator.u

    update_log!(integrator)
end

function update_log!(integrator::Integrator)
    @unpack log_config, log_book = integrator
    @unpack u, t, dt = log_config
    
    # STORE INTERMEDIATE STATES OF THE SOLUTION 
    u ? push!(log_book.u_log, copy(integrator.u)) : nothing

    # STORE INTERMEDIATE TIMES OF THE SIMULATION
    t ? push!(log_book.t_log, integrator.t) : nothing 

    # STORE INTERMEDIATE TIMESTEPS OF THE SIMULATION
    dt ? push!(log_book.dt_log, integrator.dt) : nothing

end

