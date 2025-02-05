struct OptmuGSol
    mu
    G
end

function modified_scheme_step_mu(mu::AbstractArray{T}, sol::OptForEntropySol) where T
    sol.u_approx[end] .+ mu[1:end-1]
end

function compute_modified_bounds_mu(mu::AbstractArray{T}, sol::OptForEntropySol) where T
    equation, method, domain, modifiedDataType, boundsType = sol.equation, sol.method, sol.domain, sol.modifiedDataType, sol.boundsType
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    sL, sR = get_sL(method), get_sR(method)
    u, up = sol.u_approx[end-1], sol.u_approx[end]
    # source
    z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    M_delta, m_delta = zeros(T, Nx+1), zeros(T, Nx+1)
    for j in 1:Nx
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
                .+ mu[mod1(k,Nx)]
                
        end

        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        m_delta[j+1], M_delta[j+1] = m[1], M[1]
    end
    m_delta, M_delta
end

function constraint_muG(muG::AbstractVector{T}, sol::OptForEntropySol) where T
    L = size(muG)[1]
    mu, G = muG[1:Int(L/2),:], vcat(muG[Int(L/2)+1:end,:], muG[Int(L/2)+1,:])
    equation, domain = sol.equation, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    up = sol.u_approx[end]
    cons = 0.0

    # First constraint: Negative diffusion
    up_mod = modified_scheme_step_mu(mu, sol)
    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    cons=0.0
    for j in 1:Nx
        bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        bg_delta[j] = lambd*((G[j+1] - sol.Gopt[j+1]) - (G[j] - sol.Gopt[j]))
        cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
    end
    
    # Second constraint: Consistency of numerical entropy flux
    m_delta, M_delta = compute_modified_bounds_mu(mu, sol)
    for j in 1:Nx+1
        cons += max(0, m_delta[j] - G[j])^2 + max(0, G[j] - M_delta[j])^2
    end

    cons

end


function cost_muG(muG::AbstractVector{T}) where T
    L = size(muG)[1]
    sum(muG[1:Int(L/2),:].^2)
end

function find_optimal_muG(sol::OptForEntropySol)
    #Nx, p = size(sol.u_approx[1])
    muG_init, _p = vcat(zero(sol.Gopt), sol.Gopt), 0.0
    cons(res, muG, p) = (res .= [constraint_muG(muG, sol)])
    @show constraint_muG(muG_init, sol)
    
    optprob = OptimizationFunction((x,p) -> cost_muG(x), Optimization.AutoForwardDiff(); cons = cons)
    prob = OptimizationProblem(optprob, muG_init, _p, lcons = [-1000.0], ucons = [-0.01])
    solmuG = solve(prob, Ipopt.Optimizer())
    L = size(muG_init)[1]
    OptmuGSol(solmuG.u[1:Int(L/2),:], vcat(solmuG.u[Int(L/2)+1:end,:], solmuG.u[Int(L/2)+1,:]))
end