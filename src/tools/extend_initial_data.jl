function extendInitialDataToConstant(wd::WorstData, Nx::Int)
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    u_init, z = zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))

    for k in 1:Nx
        if k<j-sL 
            u_init[k,:] = sol.worstDataMat[i,begin,:]
            z[k] = sol.worstSource[i][begin]
        elseif j-sL <= k && k <= j+sR 
            u_init[k,:] = sol.worstDataMat[i,k-(j-sL)+1,:]
            z[k] = sol.worstSource[i][k-(j-sL)+1]
        else
            u_init[k,:] = sol.worstDataMat[i,end,:]
            z[k] = sol.worstSource[i][end]
        end
    end

    u_init, z
end

function extendInitialDataToLinear(wd::WorstData, Nx::Int; boxBounds=nothing, sourceBounds=[-1.0 1.0])
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    s_init, u_init, z = zeros(Nx, p), zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))
    dx = wd.domain.dx

    if isnothing(boxBounds)
        boxBounds = zeros((p,2))
        boxBounds[:,2] = 1.0
    end

    # Reconstruction on the water height
    left_slope = (wd.worstDataMat[i,begin+1,:] .- wd.worstDataMat[i,begin,:])./dx
    right_slope = (wd.worstDataMat[i,end,:] .- wd.worstDataMat[i,end-1,:])./dx

    # # Reconstruction on the surface height
    # left_slope = (wd.worstDataMat[i,begin+1,:] .+ wd.worstSource[i][begin+1] .- wd.worstDataMat[i,begin,:] .- wd.worstSource[i][begin])./dx
    # right_slope = (wd.worstDataMat[i,end,:] .+ wd.worstSource[i][end] .- wd.worstDataMat[i,end-1,:] .- wd.worstSource[i][end-1])./dx

    # Reconstruction of the topography
    source_left_slope = (wd.worstSource[i][begin+1] .- wd.worstSource[i][begin])./dx
    source_right_slope = (wd.worstSource[i][end] .- wd.worstSource[i][end-1])./dx

    for k in 1:Nx
        if k<j-sL 
            u_init[k,:] = min.(max.(wd.worstDataMat[i,begin,:] - left_slope*dx*(j-sL-k), boxBounds[:,1]), boxBounds[:,2])
            #s_init[k,:] = min.(max.(wd.worstDataMat[i,begin,:] - left_slope*dx*(j-sL-k), boxBounds[:,1]), boxBounds[:,2])
            z[k] = min.(max.(wd.worstSource[i][begin] - source_left_slope*dx*(j-sL-k), sourceBounds[1]), sourceBounds[2])
        elseif j-sL <= k && k <= j+sR 
            u_init[k,:] = wd.worstDataMat[i,k-(j-sL)+1,:]
            z[k] = wd.worstSource[i][k-(j-sL)+1]
        else
            u_init[k,:] = min.(max.(wd.worstDataMat[i,end,:] + right_slope*dx*(k-j-sR), boxBounds[:,1]), boxBounds[:,2])
            z[k] = min.(max.(wd.worstSource[i][end] + source_right_slope*dx*(k-j-sR), sourceBounds[1]), sourceBounds[2])
        end
    end

    u_init, z
end

function extendInitialDataToK(wd::WorstData, Nx::Int)
    i = argmin(wd.worstLowDiffVec)
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)
    u_init, z = zeros(Nx, p), zeros(Nx, 1)
    j = Int(round(Nx/2))

    K = computeK(wd.modifiedDataType, extractLocalData(wd.worstDataMat[i,:,:], sL+1, sL, sR))
    Z = isnothing(wd.domain.sourceVec) ? zeros(1,p) : computeK(wd.modifiedDataType, extractLocalData(reshape(wd.worstSource[i], (sL+sR+1,1)), sL+1, sL, sR))

    for k in 1:Nx
        if j-sL ≤ k ≤ j+sR
            u_init[k,:] = wd.worstDataMat[i,k-j+sL+1,:]
            z[k,:] .= isnothing(wd.domain.sourceVec) ? 0.0 : wd.worstSource[i][k-j+sL+1]
        else
            u_init[k,:] = K
            z[k,:] = Z 
        end
    end

    z = isnothing(wd.domain.sourceVec) ? nothing : z

    u_init, z
end

function correct_extend_initial_data(wd::WorstData, args...)
    rw = wd.reducedWD
    @show rw.domain.Nx
    println(rw.eps_uz)
    #domain = rw.domain
    domain = createInterval(rw.domain)
    uz_unk = rw.uz
    p = get_unknowns_number(wd.equation)
    Nx = domain.Nx
    sL, sR = get_sL(wd.method), get_sR(wd.method)

    j = Int(round(Nx/2))
    u = zeros(Nx,p)

    # Reconstructing u from the given data
    u_unk = reshape(uz_unk[1:p*(sL+sR+1)], (sL+sR+1,p))
    K = computeK(wd.modifiedDataType, extractLocalData(u_unk, sL+1, sL, sR))
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
        Z = computeK(wd.modifiedDataType, extractLocalData(reshape(uz_unk[p*(sL+sR+1)+1:end,:], (sL+sR+1,1)), sL+1, sL, sR))
        for k in 1:Nx
            if j-sL ≤ k ≤ j+sR
                zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-j+sL+1]
            else
                zVec[k,:] = Z
            end
        end
        domain.sourceVec = zVec
    end

    # u_mid = u[j,:]
    # dx = domain.dx
    # dt = method.CFL_factor * dx / CFL_cond(equation, u) # Timestep given by CFL condition
    # z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    # # Computing mj-1/2 and Mj-1/2

    # ut = compute_u_tilde(modifiedDataType, u, j-1, sL, sR)
    # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j-1, sL, sR)
    # uh = compute_u_hat(equation.source, ut, dx, dt, j-1, domain, equation, method; zt=zt)

    # mminus, Mminus = initBounds(modifiedDataType, equation, u, j-1, sL, sR, z)
    # mminus, Mminus = updateBounds!(modifiedDataType, boundsType, equation, mminus, Mminus, ut, uh, j-1, sL, sR, Nx, dx, dt, zt)

    # if (mminus[1] > Mminus[1]) && warningsOn
    #     @warn "m-1/2 greater than M-1/2 !!!"
    # end

    # # Computing mj+1/2 and Mj+1/2

    # ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
    # uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)

    # up_mid = scheme_step(equation.source, u, dt, domain, equation, method)[j,:]

    # mplus, Mplus = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    # mplus, Mplus = updateBounds!(modifiedDataType, boundsType, equation, mplus, Mplus, ut, uh, j, sL, sR, Nx, dx, dt, zt)

    # if Mplus[1] - mplus[1] < -1e-10 && warningsOn
    #     @warn "m+1/2 greater than M+1/2 !!!"
    #     @show Mplus[1] - mplus[1]
    # end

    # @show min(get_eta(equation, u_mid; z=z)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z)[1], Mplus[1] - mplus[1])
    u, domain
end

function correct_extend_initial_data_faster(Nx::Int, wd::WorstData, args...)
    rw = wd.reducedWD
    #domain = rw.domain
    domain = createInterval(rw.domain.xmin, rw.domain.xmax, Nx, rw.domain.t0, rw.domain.Tf)
    uz_unk = rw.uz
    p = get_unknowns_number(wd.equation)
    sL, sR = get_sL(wd.method), get_sR(wd.method)

    j = Int(round(Nx/2))
    u = zeros(Nx,p)

    # Reconstructing u from the given data
    u_unk = reshape(uz_unk[1:p*(sL+sR+1)], (sL+sR+1,p))
    K = computeK(wd.modifiedDataType, extractLocalData(u_unk, sL+1, sL, sR))
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
        Z = computeK(wd.modifiedDataType, extractLocalData(reshape(uz_unk[p*(sL+sR+1)+1:end,:], (sL+sR+1,1)), sL+1, sL, sR))
        for k in 1:Nx
            if j-sL ≤ k ≤ j+sR
                zVec[k,:] .= uz_unk[p*(sL+sR+1) + k-j+sL+1]
            else
                zVec[k,:] = Z
            end
        end
        domain.sourceVec = zVec
    end

    # u_mid = u[j,:]
    # dx = domain.dx
    # dt = method.CFL_factor * dx / CFL_cond(equation, u) # Timestep given by CFL condition
    # z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    # # Computing mj-1/2 and Mj-1/2

    # ut = compute_u_tilde(modifiedDataType, u, j-1, sL, sR)
    # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j-1, sL, sR)
    # uh = compute_u_hat(equation.source, ut, dx, dt, j-1, domain, equation, method; zt=zt)

    # mminus, Mminus = initBounds(modifiedDataType, equation, u, j-1, sL, sR, z)
    # mminus, Mminus = updateBounds!(modifiedDataType, boundsType, equation, mminus, Mminus, ut, uh, j-1, sL, sR, Nx, dx, dt, zt)

    # if (mminus[1] > Mminus[1]) && warningsOn
    #     @warn "m-1/2 greater than M-1/2 !!!"
    # end

    # # Computing mj+1/2 and Mj+1/2

    # ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
    # uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)

    # up_mid = scheme_step(equation.source, u, dt, domain, equation, method)[j,:]

    # mplus, Mplus = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    # mplus, Mplus = updateBounds!(modifiedDataType, boundsType, equation, mplus, Mplus, ut, uh, j, sL, sR, Nx, dx, dt, zt)

    # if Mplus[1] - mplus[1] < -1e-10 && warningsOn
    #     @warn "m+1/2 greater than M+1/2 !!!"
    #     @show Mplus[1] - mplus[1]
    # end

    # @show min(get_eta(equation, u_mid; z=z)[1] - dt/dx *(mplus[1] - Mminus[1]) - get_eta(equation, up_mid; z=z)[1], Mplus[1] - mplus[1])
    u, domain
end