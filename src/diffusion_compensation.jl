function constraint_f_delta(f_delta::AbstractArray{T}, sol::OptForEntropySol) where T
    
    equation, method, domain, modifiedDataType, boundsType = sol.equation, sol.method, sol.domain, sol.modifiedDataType, sol.boundsType
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    sL, sR = get_sL(method), get_sR(method)
    u, up = sol.u_approx[end-1], sol.u_approx[end]
    # source
    z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))
    
    p = get_unknowns_number(equation)
    numericalFluxMat = zeros(T, (Nx+1, p))

    for i ∈ 2:Nx
        numericalFluxMat[i,:] = giveNumFlux(method, equation, u[i-1,:], u[i,:]; zL=z[i-1], zR=z[i])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, u[end,:], u[1,:]; zL=z[end], zR=z[1])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]

    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    M_delta, m_delta = zeros(T, Nx+1), zeros(T, Nx+1)
    cons=0.0
    for j in 1:Nx
        # Calcul of Bu delta
        bu_delta[j] = get_eta(equation, u[j] - lambd*(numericalFluxMat[j] - numericalFluxMat[mod1(j-1,Nx)]) - lambd*(f_delta[j+1] - f_delta[j])) - get_eta(equation, up[j])
        
        # Calcul of Bg delta
        ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        zt = isnothing(domain.sourceVec) ? zeros(T, size(ut)) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
            
                ###### This part needs to be completed for equation with source term
        Dzt = isnothing(domain.DSourceVec) ? zeros(T, size(ut)) : domain.DSourceVec
        uh = copy(ut)
        #sourceVec = isnothing(sourceTerm(method, domain, ut; z=zt)) ? zero(Nx) : sourceTerm(method, domain, ut; z=zt)
        sourceVec = sourceTerm(equation, method, domain, ut; z=zt, Dz=Dzt)
        for k in j-sL-sR+1:j+sR+sL
            uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- lambd .* (
                giveNumFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
                .-
                giveNumFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])) .+ dt * sourceVec[mod1(k, Nx),:]
                .-
                lambd*(f_delta[mod1(k, Nx)] - f_delta[mod1(k-1, Nx)])
        end

        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        m_delta[j+1], M_delta[j+1] = m[1], M[1]

        bg_delta[j] = lambd*(((M_delta[j+1]+m_delta[j+1])/2 - sol.Gopt[j+1]) - ((M_delta[j]+m_delta[j])/2 - sol.Gopt[j]))
        
        cons += max(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
        
    end
    cons
    
    #f_delta[1]
end

function f_delta_from_A(A::AbstractArray{T}, sol::OptForEntropySol) where T

    # A should be a vector
    u = sol.u_approx[end-1]
    Nx, p = size(u)
    f_delta = zeros(T, (Nx+1, p))
    for i in 1:Nx+1
        f_delta[i,:] = -0.5*A[i].* (u[mod1(i+1, Nx), :] - u[mod1(i, Nx), :])
    end
    #zero(A) .+ (1-A[1]^2)
    f_delta
end

function constraint_A(A, sol::OptForEntropySol)
    constraint_f_delta(f_delta_from_A(A, sol), sol)
end

function find_optimal_A(sol::OptForEntropySol)
    @show A_init, _p = zero(sol.Gopt).-3, 0.0
    cons(res, A, p) = (res .= [constraint_A(A, sol)])
    @show constraint_A(A_init, sol)

    #optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoForwardDiff(); cons = cons)
    optprob = OptimizationFunction((x,p) -> sum(x.^2), Optimization.AutoForwardDiff(); cons = cons)
    #optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoZygote(); cons = cons)
    #optprob = OptimizationFunction((x,p) -> abs.(x); cons = cons)
    # optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoModelingToolkit(); cons = cons)
    prob = OptimizationProblem(optprob, A_init, _p, lcons = [-Inf], ucons = [0.0])
    #solA = solve(prob, IPNewton())
    #solA = solve(prob, BFGS())
    solA = solve(prob, Ipopt.Optimizer())
end