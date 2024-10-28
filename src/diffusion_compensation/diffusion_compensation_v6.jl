# Diffusion compensation with vectorial constraints

function vec_constraint_alphaG(alphaG::AbstractArray{T}, sol::OptForEntropySol, entMethod::FVMethod) where T
    alpha, G = divideAlphaG(alphaG)
    equation, domain = sol.equation, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    u, up = sol.u_approx[end-1], sol.u_approx[end]
    cons = 0.0

    @assert sol.method.CFL_factor == entMethod.CFL_factor
    alphaMethod = MixedMethod(CFL_factor, sol.method, entMethod, alpha)

    cons = zeros(4)

    # First constraint: Negative diffusion
    up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaMethod)
    D_mod = diffusion(u, up_mod, G, dx, dt, equation, domain)
    for d in D_mod
        cons[1] += max(0, d)^2
    end
    # display(plot(sol.domain.x, D_mod))
    # plot(sol.domain.x, u)
    # display(plot!(sol.domain.x, up_mod))
    # println("After diffusion, cons="*string(cons))
    
    # Second constraint: Consistency of numerical entropy flux
    m_delta, M_delta = compute_modified_bounds_alpha(alpha, sol, alphaMethod)
    #m_delta, M_delta = compute_G_bounds(u, Nx, dx, dt, equation, domain, alphaMethod, sol.modifiedDataType, sol.boundsType)
    # plot(domain.interfaces, m_delta, label="m")
    # plot!(domain.interfaces, G, label="G")
    # plot!(domain.interfaces, M_delta, label="M")
    # display(title!("bounds of modified scheme"))
    for j in 1:Nx+1
        cons[2] += (dt/dx)^2*(max(0, m_delta[j] - G[j])^2 + max(0, G[j] - M_delta[j])^2)
    end
    # println("After consistency, cons="*string(cons))


    # Third constraint: alpha must be between specified bounds
    for a in alpha
        #cons += max(0, a-10)^2 + max(0, -9-a)^2
        cons[3] += max(0, a-1)^2 + max(0, -a)^2
    end
    # println("After alpha domain, cons="*string(cons))

    # Fourth constraint: maximum principle
    min_u, max_u = minimum(u), maximum(u)
    for j in 1:Nx
        cons[4] += max(0, min_u-up_mod[j]).^2 + max(0, up_mod[j]-max_u).^2
    end

    cons
end


function find_optimal_alphaG_vec_cons(sol::OptForEntropySol, entMethod::FVMethod, Gent)
    alphaG_init, _p = vcat(zeros(Nx).+1, Gent[1:Nx]), 0.0
    # alphaG_init, _p = vcat(zero(sol.Gopt[:end-1]).+1, Gent[:end-1]), 0.0
    alphaG_init = reshape(alphaG_init, size(alphaG_init)[1])

    cons(res, alphaG, p) = (res .= constraint_alphaG(alphaG, sol, entMethod))
    @show constraint_alphaG(alphaG_init, sol,  entMethod)

    cost(alphaG) = cost_alphaG(alphaG, sol, entMethod)

    optprob = OptimizationFunction((x,p) -> cost(x), Optimization.AutoForwardDiff(); cons = cons)
    prob = OptimizationProblem(optprob, alphaG_init, _p, lcons = [-1000.0,-1000.0,-1000.0,-1000.0], ucons = [0.0,0.0,0.0,0.0])
    solAlphaG = solve(prob, Ipopt.Optimizer(); max_iter=100, tol=1e-2)
    
    alphaSol, GSol = divideAlphaG(solAlphaG.u)
    equation, domain = sol.equation, sol.domain
    dt = sol.dt_vec[end]
    u = sol.u_approx[end-1]
    alphaMethod = MixedMethod(CFL_factor, sol.method, entMethod, alphaSol)
    up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaMethod)

    OptAlphaGSol(alphaSol, GSol, entMethod, up_mod)
end