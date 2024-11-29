# A priori quantification of numerical diffusion which could easily be extended to a multidimensional system

function compute_u_tilde_multidim(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int)
    
    K = computeK(KFun, extractLocalData_multidim(u, j, sL, sR))
    Nx= size(u)[1]
    ut = zero(u)
    
    for i in 1:Nx
        # if (j - sL <= i && i <= j + sR) || (((j - sL) % Nx <= i) && (j - sL < 1)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
        if (j - sL <= i && i <= j + sR) || ((mod1(j - sL, Nx) <= i) && (j - sL <= 0)) || ((i <= mod1(j + sR, Nx)) & (j + sR > Nx))
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

function isBetween(a::Int, b::Int, k::Int, N::Int)
    if a ≤ k && k ≤ b
        return true
    elseif a ≤ k-N ≤ b 
        return true
    elseif a ≤ k+N ≤ b
        return true
    end
    false
end

function compute_u_tilde_multidim(::AsymmetricModifiedData, u, j::Int, sL::Int, sR::Int)

    Nx = size(u)[1]
    ut = zero(u)

    ileft, iright = j - sL, j + sR
    # jleft, jright = j - sL - sR + 1, j + sL + sR - 1
    jleft, jright = j - sL - sR, j + sL + sR
    # @show j
    # @show ileft, iright
    # @show jleft, jright

    for k in 1:Nx

        # First, we fill the cells necessary for consistency
        if isBetween(ileft, iright, k, Nx)
            ut[k,:] = u[k,:]
        elseif isBetween(jleft, ileft, k, Nx)
            ut[k,:] = u[mod1(ileft, Nx),:]
            #println("ok left")
        elseif isBetween(iright, jright, k, Nx)
            ut[k,:] = u[mod1(iright, Nx),:]
            #println("ok right")
        
        # Then we extend on both sides if cell value isn't already defined (but in practice, this step could be supressed)
        elseif 1 ≤ k && k ≤ jleft
            ut[k,:] = u[mod1(ileft, Nx),:]
        elseif jright ≤ k && k ≤ Nx
            ut[k,:] = u[mod1(iright, Nx),:]
        end
    end

    # if j==1
    #     @show ileft, iright, jleft, jright
    # end

    # for k in 1:Nx
    #     #@show k

    #     if (
    #         (ileft <= k && iright >= k)
    #         || (mod1(ileft, Nx) <= k && ileft < 1)
    #         || (mod1(iright, Nx) >= k && iright > Nx)
    #     )
    #     # if (
    #     #     (mod1(ileft, Nx) <= k && mod1(iright, Nx) >= k)
    #     #     || (mod1(ileft, Nx) > mod1(iright, Nx)
    #     #     && (k<= mod1(iright, Nx) || k>= mod1(ileft,Nx)))
    #     # )
    #         ut[k,:] = u[k,:]
    #         # @show mod1(ileft, Nx) <= k && mod1(iright,Nx) >= k
    #         # @show mod1(ileft, Nx) > mod1(iright, Nx)
    #         # @show k<= iright
    #         # @show k>= ileft
    #         # @show ileft <= k && iright >= k
    #         # @show ileft % Nx <= k && ileft < 1
    #         # @show iright % Nx >= k && iright > Nx
    #         #println("Keep value")
    #     elseif (jleft <= k && ileft > k) || (
    #         mod1(jleft, Nx) <= k && jleft < 1 && k < mod1(ileft, Nx)
    #     )
    #         ut[k,:] = u[mod1(ileft, Nx),:]
    #         #println("left value")
    #     elseif (iright < k && jright >= k) || (
    #         mod1(jright, Nx) >= k && jright > Nx && k > mod1(iright, Nx)
    #     )
    #         ut[k,:] = u[mod1(iright, Nx),:]
    #         #println("rigth")
    #     elseif jright < k
    #         ut[k,:] = u[mod1(iright, Nx),:]
    #         #println("right")
    #     elseif jleft > k
    #         ut[k,:] = u[mod1(ileft, Nx),:]
    #         #println("left")
    #     end
    # end
    ut
end

function compute_u_hat_multidim(::NullSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=nothing)

    # uh = copy(ut)
    # Nx, p = size(ut)

    # sL, sR = get_sL(method), get_sR(method)

    # # numericalFluxMat = giveNumFlux(zb, method, equation, v; domain=domain)

    # for k in j-sL-sR:j+sR+sL

    #     uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
    #         numFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :])
    #         .-
    #         numFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :])
    #     )

    # end
    uh = scheme_step(equation.source, ut, dt, domain, equation, method)
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

function init_l(::SymmetricModifiedData, u, args...)
    zero(eltype(u))
end

function update_l!(::SymmetricModifiedData, equation::Equation, l, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))
    for k in j-sL-sR:j-1
        l += get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1] - get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1]
    end
        
    for k in j+1:j+sL+sR
        l += get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1] - get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1]
    end
    l
end

function init_l(::AsymmetricModifiedData, u, equation::Equation, j, sL, sR, dx, dt, z=nothing)
    Nx = size(u)[1]
    z = nothing
    zL = isnothing(z) ? nothing : z[mod1(j + sR, Nx),:]
    zR = isnothing(z) ? nothing : z[mod1(j - sL, Nx),:]
    dt/dx * (get_G(equation, u[mod1(j + sR, Nx),:], zR)[1] - get_G(equation, u[mod1(j - sL, length(u)),:], zL)[1])
end

function update_l!(::AsymmetricModifiedData, equation::Equation, l, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))
    for k in j-sL-sR+1:j-1
        l += get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1] - get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1]
    end
        
    for k in j+1:j+sL+sR-1
        l += get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1] - get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])[1]
    end
    l
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
        @test ut[j,:]==u[j,:]
        for k in j-sL:j+sR 
            @test ut[mod1(k,Nx)] == u[mod1(k,Nx)]
        end
        @test uh[j,:] == scheme_step(equation.source, ut, dt, domain, equation, method)[j,:]
        L = get_eta(equation, ut[j,:], zt[j,:]) .- get_eta(equation, uh[j,:], zt[j,:])

        # Calculating lj+1/2
        # p = get_unknowns_number(equation)
        # l = zeros((p,))
        l = init_l(modifiedDataType, u, equation::Equation, j, sL, sR, dx, dt, z)
        l = update_l!(modifiedDataType, equation, l, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        # if j==1
        #     plot(ut[:,1], label="ut")
        #     display(plot!(uh[:,1], label="uh"))
        # end

        # plot(ut[:,1], label="ut")
        # plot!(uh[:,1], label="uh")
        # display(title!("j = "*string(j)))

        # for k in j-sL-sR:j-1
        #     l .+= get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])
        #     # if j==1 || j==2
        #     #     println("New iteration")
        #     #     @show k
        #     #     @show mod1(k, Nx)
        #     #     @show l 
        #     #     @show get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
        #     # end
        # end
            
        # for k in j+1:j+sL+sR
        #     l .+= get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:])
        #     # if j==1 || j==2
        #     #     println("New iteration")
        #     #     @show k
        #     #     @show mod1(k, Nx)
        #     #     @show l 
        #     #     @show get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:])
        #     # end
        # end

        l_vec[j+1], L_vec[j+1] = l[1], L[1]

    end

    l_vec[1], L_vec[1] = l_vec[end], L_vec[end]

    l_vec, L_vec
end

function diffusion_a_priori_multidim(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK_multidim(get_sL(method), get_sR(method)), boundsType::BoundsType=NormalBounds())
    
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

    dlL = L_vec[begin+1:end].-l_vec[begin+1:end]
    display(plot(domain.x, dlL, label="L-l"))
    @show min.(dlL,0.0)

    @show sum(L_vec[begin+1:end])
    @show sum(l_vec[begin+1:end])

    # D_low = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*l_vec[begin+1:end]
    # D_up = [get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])].+ dt_vec[end]/dx*L_vec[begin+1:end]

    
    D_base = [get_eta(equation, u_approx[end][i,:], z[i])[1] - get_eta(equation, u_approx[end-1][i,:], z[i])[1] for i in 1:length(u_approx[end-1][:,1])]
    
    @test sum(L_vec[begin+1:end]) == -sum(D_base)

    D_low = D_base .+ l_vec[begin+1:end]
    D_up = D_base .+ L_vec[begin+1:end]


    # println("Normalisation")
    # @show sum([get_eta(equation, u_approx[end][i,:]; z=z[i])[1] - get_eta(equation, u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(u_approx[end-1][:,1])])
    # @show sum(D_low)

    # The normalisation is chosen to keep the same amount of diffusion up to a sign 
    # We want to keep the same sign as the pre-normalised quantity

    D_low_norm = D_low * sum(D_base) / sum(D_low)
    D_up_norm = D_up * sum(D_base) / abs(sum(D_up))

    # # Test d'une autre normalisation (par fenêtre glissante)
    # J = Int(Nx/10)
    # D_low_norm=zero(D_low)
    # for j in 1:Nx
    #     s1, s2 = 0.0, 0.0
    #     for k in j-J+1:j+J
    #         s1 += (D_low[mod1(k,Nx)]+D_up[mod1(k,Nx)])/2
    #         s2 += D_low[mod1(k,Nx)]
    #     end
    #     D_low_norm[j] = D_low[j]*s1/s2
    # end

    PrioriDiffSol(domain, equation, method, modifiedDataType, boundsType, u_approx, dt_vec, l_vec, L_vec, D_low, D_up, D_low_norm, D_up_norm, compute_D_CL(D_low, D_up, D_base)...)
end