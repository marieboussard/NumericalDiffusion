struct LsGsol
    domain::Domain
    equation::Equation
    scheme::FVScheme
    u0
    Gc::NumEntFlux
    modifiedDataType::ModifiedDataType
    l
    L
    is_solved_and_feasible
    Gopt
    cost
end

function compute_delta_G(gamma)
    Nx = length(gamma)
    delta_G = zero(gamma)
    for j in 1:Nx
        delta_G = gamma[mod1(j, Nx)] - gamma[mod1(j - 1, Nx)]
    end
    delta_G
end

function lsEntropicConsistentG(domain::Domain, equation::Equation, scheme::FVScheme, u0, Gc::NumEntFlux; modifiedDataType::ModifiedDataType=meanK_multidim(get_sL(scheme), get_sR(scheme)))
    #dt = computeTimeStep(method, equation, domain, u0)
    dt = domain.Tf
    @show l, L = compute_multidim_bounds(u0, domain.Nx, domain.dx, dt, equation, domain, scheme, modifiedDataType)
    @show gc = vecNumFlux(equation.source, Gc, equation, u0; dt=dt, domain=domain)[end-1]
    gamma_init, p0 = zeros(domain.Nx), 0.0
    # cons(res, x, p) = (res .= compute_delta_G(x))
    # prob = Optimization.OptimizationProblem((x,p)->norm(x.-gc, 2), gamma_init, p0, l[begin+1:end], L[begin+1:end])
    # sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 1000,
    # maxtime = 1000.0)
    # optprob = OptimizationFunction((x,p)->norm(x.-gc, 2), Optimization.AutoForwardDiff(); cons = cons)
    # prob = OptimizationProblem(optprob, gamma_init, p0, lcons = l[begin+1:end], ucons = L[begin+1:end])
    # solAG = solve(prob, Ipopt.Optimizer())
    # optprob = OptimizationFunction((x,p)->norm(x.-gc, 2), Optimization.AutoForwardDiff(), cons = cons)
    # prob = OptimizationProblem(optprob, gamma_init, p0, lcons = l[begin+1:end], ucons = L[begin+1:end])
    # sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())
    # sol

    # Resolution with OSQP
    lamb = dt / domain.dx
    model = Model(OSQP.Optimizer)
    @variable(model, gamma[1:domain.Nx])
    for j in 1:domain.Nx
        @constraint(model, l[j+1] ≤ lamb * (gamma[j] - gamma[mod1(j - 1, domain.Nx)]))
        @constraint(model, lamb * (gamma[j] - gamma[mod1(j - 1, domain.Nx)]) ≤ L[j+1])
    end
    @objective(model, Min, sum((gamma .- gc) .^ 2))

    #print(model)
    optimize!(model)
    Gopt = zeros(typeof(value(gamma[1])), domain.Nx + 1)
    Gopt[begin+1:end] = value.(gamma)
    Gopt[begin] = Gopt[end]
    LsGsol(domain, equation, scheme, u0, Gc, modifiedDataType, l, L, is_solved_and_feasible, Gopt, objective_value(model))
end