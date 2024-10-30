# A priori quantification of numerical diffusion which could easily be extended to a multidimensional system

function compute_u_tilde_multidim(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int)
    
    K = computeK(KFun, extractLocalData(u, j, sL, sR))
    Nx= size(u)[1]
    ut = zero(u)
    
    for i in 1:Nx
        # if (j - sL <= i && i <= j + sR) || (((j - sL) % Nx <= i) && (j - sL < 1)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
        if (j - sL <= i && i <= j + sR) || (((j - sL) % Nx <= i) && (j - sL < 0)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
        ut[i,:] = u[i,:]
            # if j==1
            #     println(i)
            #     @show (j - sL <= i && i <= j + sR)
            #     @show (((j - sL) % Nx <= i) && (j - sL < 1))
            #     @show ((i <= (j + sR) % Nx) & (j + sR > Nx))
            # end
        else
            ut[i,:] .= K
        end
    end

    ut

end

function compute_u_hat_multidim(::NullSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=nothing)

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    # numericalFluxMat = giveNumFlux(zb, method, equation, v; domain=domain)

    for k in j-sL-sR:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            numFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :])
            .-
            numFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :])
        )

    end
    #scheme_step(ns, ut, dt, domain, equation, method)
    return uh
end

function compute_u_hat_multidim(::ZbSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=zero(ut))

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    sourceVec = sourceTerm(equation, method, domain, ut; z=zt)

    for k in j-sL-sR:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            giveNumFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            .-
            giveNumFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])) .+ dt * sourceVec[mod1(k, Nx),:]
        

    end

    return uh
end


function compute_multidim_bounds(u, Nx, dx, dt, equation::Equation, domain::Domain, method::FVMethod, modifiedDataType::ModifiedDataType)
    l_vec, L_vec = zeros(Nx + 1), zeros(Nx + 1)
    sL, sR = get_sL(method), get_sR(method)

    # source
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    for j in 1:Nx
        ut = compute_u_tilde_multidim(modifiedDataType, u, j, sL, sR)
        zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde_multidim(modifiedDataType, z, j, sL, sR)
        uh = compute_u_hat_multidim(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)

        # Calculating Lj+1/2
        L = get_eta(equation, ut[j,:]; z=zt[j,:]) .- get_eta(equation, uh[j,:]; z=zt[j,:])

        # Calculating lj+1/2
        p = get_unknowns_number(equation)
        l = zeros((p,))

        if j==1
            plot(ut[:,1], label="ut")
            display(plot!(uh[:,1], label="uh"))
        end

        for k in j-sL-sR:j-1
            l .+= get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
            # if j==1 || j==2
            #     println("New iteration")
            #     @show k
            #     @show mod1(k, Nx)
            #     @show l 
            #     @show get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
            # end
        end
            
        for k in j+1:j+sL+sR
            l .+= get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
            # if j==1 || j==2
            #     println("New iteration")
            #     @show k
            #     @show mod1(k, Nx)
            #     @show l 
            #     @show get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
            # end
        end

        l_vec[j+1], L_vec[j+1] = l[1], L[1]

    end

    l_vec[1], L_vec[1] = l_vec[end], L_vec[end]

    l_vec, L_vec
end

function diffusion_a_priori_multidim(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(get_sL(method), get_sR(method)), boundsType::BoundsType=NormalBounds())
    
    Nx, dx = domain.Nx, domain.dx

    FVsol = fv_solve(domain, u_init, equation, method)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    # source
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    l_vec, L_vec = compute_multidim_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, method, modifiedDataType)

    solEnt = optimize_for_entropy(u_init, domain, equation, method)
    plot(domain.x, l_vec[begin+1:end], label="l")
    display(plot!(domain.x, (solEnt.Gopt[begin+1:end].-solEnt.Gopt[begin:end-1])/domain.dx*dt_vec[end], label="dG"))
    display(plot!(domain.x, L_vec[begin+1:end], label="L"))

    # D_low = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*l_vec[begin+1:end]
    # D_up = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*L_vec[begin+1:end]

    D_low = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ l_vec[begin+1:end]
    D_up = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ L_vec[begin+1:end]


    # println("Normalisation")
    # @show sum([get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])])
    # @show sum(D_low)

    # The normalisation is chosen to keep the same amount of diffusion up to a sign 
    # We want to keep the same sign as the pre-normalised quantity

    #D_low_norm = D_low * sum([get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / abs(sum(D_low))
    D_up_norm = D_up * sum([get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])]) / abs(sum(D_up))

    # Test d'une autre normalisation (par fenêtre glissante)
    J = Int(Nx/10)
    D_low_norm=zero(D_low)
    for j in 1:Nx
        s1, s2 = 0.0, 0.0
        for k in j-J+1:j+J
            s1 += (D_low[mod1(k,Nx)]+D_up[mod1(k,Nx)])/2
            s2 += D_low[mod1(k,Nx)]
        end
        D_low_norm[j] = D_low[j]*s1/s2
    end

    PrioriDiffSol(domain, equation, method, modifiedDataType, boundsType, u_approx, dt_vec, l_vec, L_vec, D_low, D_up, D_low_norm, D_up_norm)
end