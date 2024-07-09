function epsilon(u_unk, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType)

    p = get_unknowns_number(equation)
    Nx = domain.Nx

    u = zeros(Nx,p)

    #u = reshape(u, (Nx,p))
    sL, sR = get_sL(method), get_sR(method)
    j = sL + sR + 1

    for k in sR+1:sL+2*sR+1
        u[k,:] .= u_unk[k-sR]
    end

    u_mid = u[j]
    dx = domain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u) # Timestep given by CFL condition

    # Computing mj-1/2 and Mj-1/2

    ut = compute_u_tilde(modifiedDataType, u, j-1, sL, sR)
    uh = compute_u_hat(ut, dx, dt, j-1, domain, equation, method)

    # source
    z = zeros(domain.Nx,1)
    for i in eachindex(domain.x) z[i]=zb(equation.source, domain.x[i]) end
    zt = compute_z_tilde(equation.source, modifiedDataType, domain, j-1, sL, sR)

    mminus, Mminus = initBounds(modifiedDataType, equation, u, j-1, sL, sR, z)
    mminus, Mminus = updateBounds!(modifiedDataType, boundsType, equation, mminus, Mminus, ut, uh, j-1, sL, sR, Nx, dx, dt, zt)

    #@show mminus, Mminus
    if mminus > Mminus
        @warn "m-1/2 greater than M-1/2 !!!"
    end

    # Computing mj+1/2 and Mj+1/2

    ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    uh = compute_u_hat(ut, dx, dt, j, domain, equation, method)

    #@show up_mid = uh[j]
    up_mid = scheme_step(equation.source, u, dt, domain, equation, method)[j]

    # source
    z = zeros(domain.Nx,1)
    for i in eachindex(domain.x) z[i]=zb(equation.source, domain.x[i]) end
    zt = compute_z_tilde(equation.source, modifiedDataType, domain, j, sL, sR)

    mplus, Mplus = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    mplus, Mplus = updateBounds!(modifiedDataType, boundsType, equation, mplus, Mplus, ut, uh, j, sL, sR, Nx, dx, dt, zt)

    #@show mplus, Mplus
    if Mplus[1] - mplus[1] < -1e-10
        @warn "m+1/2 greater than M+1/2 !!!"
        @show Mplus[1] - mplus[1]
    end

    # @show get_eta(equation, u_mid; z=z)
    # @show get_eta(equation, up_mid; z=z)[1]
    # @show mplus[1] - Mminus[1]

    res = min(get_eta(equation, u_mid; z=z)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z)[1], Mplus[1] - mplus[1])

    #@show res

    res

end

function find_worst_initial_data(Nx, u_init, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType=NormalBounds())
    
    sL, sR = get_sL(method), get_sR(method)
    domain = createInterval(-1.0/Nx, 1.0/Nx, 2*sL+2*sR+1, 0.0, 1.0)
    p = get_unknowns_number(equation)
    
    #u_init = zeros(sL+sR+1, p) .+1
    #u_init = zeros((sL+sR+1)*p) .+ 2

    lower = -1 .+ u_init
    upper = 1 .+ u_init

    #df = TwiceDifferentiable(u -> epsilon(u, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType), u_init)

    #dfc = TwiceDifferentiableConstraints(lower, upper)

    #sol = optimize((u,p) -> epsilon(u, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType), dfc, u_init; iterations=10000)
    
    prob = OptimizationProblem((u, p) -> epsilon(u, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType), u_init, p=u_init, lb = lower, ub = upper)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())

    worstData, worstLowDiff = sol.minimizer, sol.minimum
end


function iterate_WID(Nx, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType=NormalBounds(), nb_it=100)

    box_size = 10.0
    p = get_unknowns_number(equation)
    sL, sR = get_sL(method), get_sR(method)
    res = zeros(nb_it)
    for k in 1:nb_it
        u_init = zeros((sL+sR+1)*p) .+ rand()*box_size
        _, res[k] = find_worst_initial_data(Nx, u_init, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType)
    end
    res
end

