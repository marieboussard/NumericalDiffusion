struct ReducedWD
    uz
    eps_uz
    domain::Domain
end
struct WorstData
    initDataMat
    initSource
    worstDataMat
    worstLowDiffVec
    worstSource
    equation::Equation
    method::FVMethod
    domain::Domain
    modifiedDataType::ModifiedDataType
    reducedWD
end

function plotWorstWD(wd::WorstData)
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    for k in 1:p
        scatter(wd.initDataMat[i,:,k], label="Initial")
        scatter!(wd.worstDataMat[i,:,k], label="Worst")
        xlabel!("x")
        title!("Worst Initial Data for "*get_name(wd.equation)*" with "*get_name(wd.method))
        display(ylabel!("Unknow number "*string(k)))
    end
    if !isnothing(wd.worstSource[i])
        scatter(wd.initSource[i], label="Initial")
        scatter!(wd.worstSource[i], label="Worst")
        xlabel!("x")
        ylabel!("Source")
        display(title!("Worst Initial Source for "*get_name(wd.equation)*" with "*get_name(wd.method)))
    end
    println("Worst value found for epsilon : "*string(wd.worstLowDiffVec[i]))
end

function plotWorstWD(wd::WorstData, ::SaintVenant, plotMode::PlottingMode=DisplayMode())
    i = argmin(wd.worstLowDiffVec)

    pltA = []

    plt = createPlot(plotMode)

    plot!(wd.initDataMat[i,:,1] .+ wd.initSource[i], label="Initial", markershape = :circ, mc=:blue3, lc=:blue3, lw=2, ms=8)
    plot!(wd.worstDataMat[i,:,1] .+ wd.worstSource[i], label="Worst", markershape = :circ, mc=:red3, lc=:red3, lw=2, ms=8)
    plot!(wd.initSource[i], label="Initial topography", markershape = :x, markerstrokewidth=2, mc=:blue3, lc=:blue3, lw=2, ms=8)
    plot!(wd.worstSource[i], label="Worst topography", markershape = :x, markerstrokewidth=2, mc=:red3, lc=:red3, lw=2, ms=8)
    xlabel!("x")
    title!("Worst Initial Data for Saint-Venant with "*get_name(wd.method))
    ylabel!("Water Surface")

    push!(pltA, plt)

    plt2 = createPlot(plotMode)

    plot!(wd.initDataMat[i,:,2], label="Initial", markershape = :circ, mc=:blue3, lc=:blue3, lw=2, ms=8)
    plot!(wd.worstDataMat[i,:,2], label="Worst", markershape = :circ, mc=:red3, lc=:red3, lw=2, ms=8)
    xlabel!("x")
    title!("Worst Initial Data for Saint-Venant with "*get_name(wd.method))
    ylabel!("Water Flow")

    push!(pltA, plt2)

    display(assemblePlot(plotMode, pltA))

    println("Worst value found for epsilon : "*string(wd.worstLowDiffVec[i]))
end


function epsilon(uz_unk, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType, warningsOn=false)

    p = get_unknowns_number(equation)
    Nx = domain.Nx
    sL, sR = get_sL(method), get_sR(method)
    #j = sL + sR + 1
    j = Int(round(Nx/2))

    newDomain = createInterval(domain.xmin, domain.xmax, domain.Nx, domain.t0, domain.Tf)

    
    # Reconstructing u from the given data
    u = zeros(Nx,p)
    u_unk = reshape(uz_unk[1:p*(sL+sR+1)], (sL+sR+1,p))
    # for k in sR+1:sL+2*sR+1
    #     u[k,:] .= u_unk[k-sR,:]
    # end
    K = computeK(modifiedDataType, extractLocalData(u_unk, sL+1, sL, sR))
    # for k in 1:Nx
    #     if sR+2 ≤ k ≤ sL+2*sR+2
    #         u[k,:] .= u_unk[k-sR-1,:]
    #     else
    #         u[k,:] = K
    #     end
    # end
    for k in 1:Nx
        if j-sL ≤ k ≤ j+sR
            u[k,:] = u_unk[k-j+sL+1,:]
        else
            u[k,:] = K
        end
    end

    # Reconstructing the source from the given data
    if length(uz_unk) > p*(sL+sR+1)
        zVec = zeros(Nx,1)
        # for k in sR+1:sL+2*sR+1
        #     zVec[k] = uz_unk[p*(sL+sR+1) + k-sR]
        # end
        Z = computeK(modifiedDataType, extractLocalData(reshape(uz_unk[p*(sL+sR+1)+1:end,:], (sL+sR+1,1)), sL+1, sL, sR))
        # for k in 1:Nx
        #     if sR+2 ≤ k ≤ sL+2*sR+2
        #         zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-sR-1]
        #     else
        #         zVec[k,:] = Z
        #     end
        # end
        for k in 1:Nx
            if j-sL ≤ k ≤ j+sR
                zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-j+sL+1]
            else
                zVec[k,:] = Z
            end
        end
        newDomain.sourceVec = zVec
    end

    u_mid = u[j,:]
    dx = newDomain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u) # Timestep given by CFL condition
    z = isnothing(newDomain.sourceVec) ? zeros((Nx, 1)) : reshape(newDomain.sourceVec, (Nx,1))

    # Computing mj-1/2 and Mj-1/2

    ut = compute_u_tilde(modifiedDataType, u, j-1, sL, sR)
    zt = isnothing(newDomain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j-1, sL, sR)
    uh = compute_u_hat(equation.source, ut, dx, dt, j-1, newDomain, equation, method; zt=zt)

    mminus, Mminus = initBounds(modifiedDataType, equation, u, j-1, sL, sR, z)
    mminus, Mminus = updateBounds!(modifiedDataType, boundsType, equation, mminus, Mminus, ut, uh, j-1, sL, sR, Nx, dx, dt, zt)

    if (mminus[1] > Mminus[1]) && warningsOn
        @warn "m-1/2 greater than M-1/2 !!!"
    end

    # Computing mj+1/2 and Mj+1/2

    ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    zt = isnothing(newDomain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
    uh = compute_u_hat(equation.source, ut, dx, dt, j, newDomain, equation, method; zt=zt)

    up_mid = scheme_step(equation.source, u, dt, newDomain, equation, method)[j,:]

    mplus, Mplus = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    mplus, Mplus = updateBounds!(modifiedDataType, boundsType, equation, mplus, Mplus, ut, uh, j, sL, sR, Nx, dx, dt, zt)

    if Mplus[1] - mplus[1] < -1e-10 && warningsOn
        @warn "m+1/2 greater than M+1/2 !!!"
        @show Mplus[1] - mplus[1]
    end
    z_mid = z[j]
    # @show uz_unk
    # @show min(get_eta(equation, u_mid; z=z)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z)[1], Mplus[1] - mplus[1])

    min(get_eta(equation, u_mid; z=z_mid)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z_mid)[1], Mplus[1] - mplus[1])

end

function epsilon_faster(uz_unk, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType, warningsOn=false)

    p = get_unknowns_number(equation)
    Nx = domain.Nx
    sL, sR = get_sL(method), get_sR(method)
    j = sL + sR + 1

    newDomain = createInterval(domain.xmin, domain.xmax, domain.Nx, domain.t0, domain.Tf)

    
    # Reconstructing u from the given data
    u = zeros(Nx,p)
    u_unk = reshape(uz_unk[1:p*(sL+sR+1)], (sL+sR+1,p))
    # for k in sR+1:sL+2*sR+1
    #     u[k,:] .= u_unk[k-sR,:]
    # end
    K = computeK(modifiedDataType, extractLocalData(u_unk, sL+1, sL, sR))
    # for k in 1:Nx
    #     if sR+2 ≤ k ≤ sL+2*sR+2
    #         u[k,:] .= u_unk[k-sR-1,:]
    #     else
    #         u[k,:] = K
    #     end
    # end
    for k in 1:Nx
        if j-sL ≤ k ≤ j+sR
            u[k,:] = u_unk[k-j+sL+1,:]
        else
            u[k,:] = K
        end
    end

    # Reconstructing the source from the given data
    if length(uz_unk) > p*(sL+sR+1)
        zVec = zeros(Nx,1)
        # for k in sR+1:sL+2*sR+1
        #     zVec[k] = uz_unk[p*(sL+sR+1) + k-sR]
        # end
        Z = computeK(modifiedDataType, extractLocalData(reshape(uz_unk[p*(sL+sR+1)+1:end,:], (sL+sR+1,1)), sL+1, sL, sR))
        # for k in 1:Nx
        #     if sR+2 ≤ k ≤ sL+2*sR+2
        #         zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-sR-1]
        #     else
        #         zVec[k,:] = Z
        #     end
        # end
        for k in 1:Nx
            if j-sL ≤ k ≤ j+sR
                zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-j+sL+1]
            else
                zVec[k,:] = Z
            end
        end
        newDomain.sourceVec = zVec
    end

    u_mid = u[j,:]
    dx = newDomain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u) # Timestep given by CFL condition
    z = isnothing(newDomain.sourceVec) ? zeros((Nx, 1)) : reshape(newDomain.sourceVec, (Nx,1))

    # Computing mj-1/2 and Mj-1/2

    ut = compute_u_tilde(modifiedDataType, u, j-1, sL, sR)
    zt = isnothing(newDomain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j-1, sL, sR)
    uh = compute_u_hat(equation.source, ut, dx, dt, j-1, newDomain, equation, method; zt=zt)

    mminus, Mminus = initBounds(modifiedDataType, equation, u, j-1, sL, sR, z)
    mminus, Mminus = updateBounds!(modifiedDataType, boundsType, equation, mminus, Mminus, ut, uh, j-1, sL, sR, Nx, dx, dt, zt)

    if (mminus[1] > Mminus[1]) && warningsOn
        @warn "m-1/2 greater than M-1/2 !!!"
    end

    # Computing mj+1/2 and Mj+1/2

    ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    zt = isnothing(newDomain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
    uh = compute_u_hat(equation.source, ut, dx, dt, j, newDomain, equation, method; zt=zt)

    up_mid = scheme_step(equation.source, u, dt, newDomain, equation, method)[j,:]

    mplus, Mplus = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    mplus, Mplus = updateBounds!(modifiedDataType, boundsType, equation, mplus, Mplus, ut, uh, j, sL, sR, Nx, dx, dt, zt)

    if Mplus[1] - mplus[1] < -1e-10 && warningsOn
        @warn "m+1/2 greater than M+1/2 !!!"
        @show Mplus[1] - mplus[1]
    end
    z_mid = z[j]
    # @show uz_unk
    # @show min(get_eta(equation, u_mid; z=z)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z)[1], Mplus[1] - mplus[1])

    min(get_eta(equation, u_mid; z=z_mid)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z_mid)[1], Mplus[1] - mplus[1])

end

function find_worst_initial_data(u_init, lower, upper, equation::Equation, method::FVMethod, domain::Domain; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType=NormalBounds())
    prob = OptimizationProblem((u, p) -> epsilon_faster(u, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType), u_init, p=u_init, lb = lower, ub = upper)
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(), abstol=1e-8)
    sol.minimizer, sol.minimum

    # sol = optimize(u -> epsilon(u, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType), u_init)
    # sol.minimizer, sol.minimum
end

function initialize_WID(::NullSource, u_init, u_low, u_up, args...)
    u_init, u_low, u_up
end
function initialize_WID(::ZbSource, u_init, u_low, u_up, sL, sR, sourceBounds)
    z_init = zeros(sL+sR+1)
    z_low, z_up = zeros(sL+sR+1) .+ sourceBounds[1], zeros(sL+sR+1) .+ sourceBounds[2]
    for i in eachindex(z_init)
        z_init[i] += sourceBounds[1] + rand()*(sourceBounds[2] - sourceBounds[1])
    end
    vcat(u_init, z_init), vcat(u_low, z_low), vcat(u_up, z_up)
end

function iterate_WID(xmin, xmax, Nx, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1,1), boundsType::BoundsType=NormalBounds(), nb_it::Int=10, boxBounds=nothing, sourceBounds=[-1.0 1.0])

    dx = (xmax - xmin) / Nx
    p = get_unknowns_number(equation)
    sL, sR = get_sL(method), get_sR(method)
    #domain = createInterval(-1.0/Nx, 1.0/Nx, 2*sL+2*sR+1, 0.0, 1.0)
    #####domain = createInterval(xmin, xmax, Nx, 0.0, 1.0)
    domain = createInterval(0.0, dx*(2*sL+2*sR+3), 2*sL+2*sR+3, 0.0, 1.0)
    #domain = createInterval(0.0, dx*(2*sL+2*sR+3), 2*sL+2*sR+3, 0.0, 0.1)
    #domain = createInterval(xmin, xmax, Nx, 0.0, 1.0)
    initDataMat, worstDataMat, worstLowDiffVec = zeros((nb_it, sL+sR+1, p)), zeros((nb_it, sL+sR+1, p)), zeros(nb_it)
    initSource, worstSource = [], []

    if isnothing(boxBounds)
        boxBounds = zeros((p,2))
        boxBounds[:,2] = 1.0
    end

    rw=nothing

    for k in 1:nb_it

        domain.sourceVec = nothing

        # Generating a random componentwise constant initial data inside the specified bounds
        u_init = zeros((sL+sR+1, p))
        u_low, u_up = zero(u_init), zero(u_init)
        for i in 1:p
            for k in 1:sL+sR+1
                u_init[k,i] = boxBounds[i,1] + (boxBounds[i,2] - boxBounds[i,1])*rand()
                u_low[k,i], u_up[k,i] = boxBounds[i,1], boxBounds[i,2]
            end
        end
        initDataMat[k,:,:] = u_init
        
        u_init = reshape(u_init, (sL+sR+1)*p)
        u_low, u_up = reshape(u_low, (sL+sR+1)*p), reshape(u_up, (sL+sR+1)*p)

        # Generating a random initial source parameter inside the specified bounds
        uz_init, uz_low, uz_up = initialize_WID(equation.source, u_init, u_low, u_up, sL, sR, sourceBounds)
        length(uz_init) > (sL+sR+1)*p ? push!(initSource, uz_init[(sL+sR+1)*p+1:end]) : push!(initSource, nothing)

        worstData, worstLowDiffVec[k] = find_worst_initial_data(uz_init, uz_low, uz_up, equation, method, domain; modifiedDataType=modifiedDataType, boundsType=boundsType)
        
        # @show worstData
        # @show domain

        #@show epsilon(worstData, domain, equation, method; modifiedDataType=modifiedDataType, boundsType=boundsType)
        
        worstDataMat[k,:,:] = reshape(worstData[1:(sL+sR+1)*p], (sL+sR+1, p))
        length(worstData) > (sL+sR+1)*p ? push!(worstSource, worstData[(sL+sR+1)*p+1:end]) : push!(worstSource, nothing)
    
        rw = ReducedWD(worstData, worstLowDiffVec, domain)
    
    end

    WorstData(initDataMat, initSource, worstDataMat, worstLowDiffVec, worstSource, equation, method, domain, modifiedDataType, rw)

end

