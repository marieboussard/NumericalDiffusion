struct CounterExample_WID
    nb_it::Int
    params
    eps
end

function counter_ex_epsilon(params, Nx, CFL_factor, equation, scheme, modifiedDataType, boundsType)

    a1, b1, a2, b2 = params

    testcase = PiecewiseLinear(a1, b1, a2, b2)
    xmin, xmax = spaceBounds(testcase)
    #### NOTE : WE DO NOT ENSURE HERE THAT ONE OF THE INTERFACES IS IN X=0
    #Nx = max(Nx, integerNx(Nx, testcase))
    domain, u0 = createOneTimestepInterval(Nx, 0.0, xmin, xmax, equation, testcase, CFL_factor)
    dx = domain.dx
    dt = domain.Tf

    m_vec, M_vec = compute_G_bounds(u0, Nx, dx, dt, equation, domain, scheme, modifiedDataType,boundsType)
    l_vec, L_vec = compute_multidim_bounds(u0, Nx, dx, dt, equation, domain, scheme, modifiedDataType)

    @show eps_1 = minimum(M_vec .- m_vec)
    @show eps_2 = maximum(l_vec .- L_vec)
    @show eps_3 = -sum(L_vec[begin+1:end])
    @show eps_4 = sum(l_vec[begin+1:end])
    max(eps_1, eps_2, eps_3, eps_4)
end


function find_counter_example(Nx, CFL_factor, equation::Equation, scheme::FVScheme; modifiedDataType::ModifiedDataType=AsymmetricModifiedData(), boundsType::BoundsType=NormalBounds(), nb_it::Int=1, paramBound=10)

    minimizers_vec = zeros(nb_it, 4)
    minima_vec = zeros(nb_it)

    for k in 1:nb_it

    # Generating a set of random parameters

        @show init_params = [rand()*paramBound for it in 1:4]
        lower = [0.0 for it in 1:4]
        upper = [paramBound for it in 1:4]

        # Minimizing the cost epsilon

        prob = OptimizationProblem((x, p) -> counter_ex_epsilon(x, Nx, CFL_factor, equation, scheme, modifiedDataType, boundsType), init_params, p=0.0, lb = lower, ub = upper)

        println("Starting resolution...")
        @time sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxtime=60, maxiters=100)#, abstol=1e-8
        minimizers_vec[k,:], minima_vec[k] = sol.minimizer, sol.minimum
        println("Resolution ended!")

    end

    i = argmin(minima_vec)

    CounterExample_WID(nb_it, minimizers_vec[i,:], minima_vec[i])

end