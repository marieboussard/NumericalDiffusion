struct OptAlphaGSol
    alpha
    G
    entMethod::FVMethod
    up
end

function divideAlphaG(alphaG::AbstractArray{T}) where T
    L = size(alphaG)[1]
    vcat(alphaG[1:Int(L/2),:], alphaG[1,]), vcat(alphaG[Int(L/2)+1:end,:], alphaG[Int(L/2)+1,:])
end

function constraint_alphaG(alphaG::AbstractArray{T}, sol::OptForEntropySol, entMethod::FVMethod) where T
    alpha, G = divide(alphaG)
    equation, domain = sol.equation, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    u, up = sol.u_approx[end-1], sol.u_approx[end]
    cons = 0.0

    @assert sol.method.CFL_factor == entMethod.CFL_factor
    alphaMethod = MixedMethod(CFL_factor, sol.method, entMethod, alpha)

    # First constraint: Negative diffusion
    up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaMethod)
    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    cons=0.0
    for j in 1:Nx
        bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        bg_delta[j] = lambd*((G[j+1] - sol.Gopt[j+1]) - (G[j] - sol.Gopt[j]))
        cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
    end
    
    # Second constraint: Consistency of numerical entropy flux
    #m_delta, M_delta = compute_modified_bounds_f_delta(f_delta, sol)
    m_delta, M_delta = compute_G_bounds(u, Nx, dx, dt, equation, domain, alphaMethod, sol.modifiedDataType, sol.boundsType)
    for j in 1:Nx+1
        cons += max(0, m_delta[j] - G[j])^2 + max(0, G[j] - M_delta[j])^2
    end

    cons
end

function cost_alphaG(alphaG::AbstractArray{T}) where T
    alpha, G = divide(alphaG)
    sum(alpha[1:end-1].^2)
end

# INPUT
# Gent is an array, the numerical entropy flux associated with the entropic method

function find_optimal_alphaG(sol::OptForEntropySol, entMethod::FVMethod, Gent)
    alphaG_init, _p = vcat(zero(sol.Gopt[:end-1]).+1, Gent[:end-1]), 0.0
    alphaG_init = reshape(alphaG_init, size(alphaG_init)[1])

    cons(res, alphaG, p) = (res .= [constraint_alphaG(alphaG, sol, entMethod)])
    @show constraint_alphaG(alphaG_init, sol,  entMethod)

    optprob = OptimizationFunction((x,p) -> cost(x), Optimization.AutoForwardDiff(); cons = cons)
    prob = OptimizationProblem(optprob, alphaG_init, _p, lcons = [-1000.0], ucons = [0.0])
    solAlphaG = solve(prob, Ipopt.Optimizer())
    
    alphaSol, GSol = divideAlphaG(solAlphaG.u)
    equation, domain = sol.equation, sol.domain
    dt = sol.dt_vec[end]
    u = sol.u_approx[end-1]
    alphaMethod = MixedMethod(CFL_factor, sol.method, entMethod, alpha)
    up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaMethod)

    OptAlphaGSol(alphaSol, GSol, entMethod, up_mod)
end