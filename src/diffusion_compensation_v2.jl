struct OptAGSol
    A
    G
end

function constraint_AG(AG::AbstractVector{T}, sol::OptForEntropySol) where T
    L = size(AG)[1]
    A, G = AG[1:Int(L/2),:], AG[Int(L/2)+1:end,:]
    equation, domain = sol.equation, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    up = sol.u_approx[end]
    cons = 0.0

    # First constraint: Negative diffusion
    f_delta = f_delta_from_A(A, sol)
    up_mod = modified_scheme_step_f_delta(A, sol)
    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    cons=0.0
    for j in 1:Nx
        bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        bg_delta[j] = lambd*((G[j+1] - sol.Gopt[j+1]) - (G[j] - sol.Gopt[j]))
        cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
    end
    
    # Second constraint: Consistency of numerical entropy flux
    m_delta, M_delta = compute_modified_bounds_f_delta(f_delta, sol)
    for j in 1:Nx+1
        cons += max(0, m_delta[j] - G[j])^2 + max(0, G[j] - M_delta[j])^2
    end

    cons

end

function cost(AG::AbstractVector{T}) where T
    L = size(AG)[1]
    sum(AG[1:Int(L/2),:].^2)
end

function find_optimal_AG(sol::OptForEntropySol)
    #Nx, p = size(sol.u_approx[1])
    @show AG_init, _p = zero(solEnt.Gopt), 0.0
    cons(res, AG, p) = (res .= [constraint_AG(AG, sol)])
    @show constraint_AG(AG_init, sol)
    
    optprob = OptimizationFunction((x,p) -> cost(x), Optimization.AutoForwardDiff(); cons = cons)
    prob = OptimizationProblem(optprob, AG_init, _p, lcons = [-1000.0], ucons = [0.0])
    solAG = solve(prob, Ipopt.Optimizer())
    L,- = size(AG)
    OptAGSol(solAG.u[1:Int(L/2),:], solAG.u[Int(L/2)+1:end,:])
end