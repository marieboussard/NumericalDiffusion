struct OptAlphaGSol{T<:AbstractFloat}
    alpha::Array{T}
    G::Array{T}
    entScheme::FVScheme
    up::Array{T}
    OptAlphaGSol{T}() where {T<:AbstractFloat} = new()
    OptAlphaGSol{T}(alpha::Array{T}, G::Array{T}, entScheme::FVScheme, up::Array{T}) where {T<:AbstractFloat} = new(alpha, G, entScheme, up)
end

OptAlphaGSol(alpha::Array{T}, G::Array{T}, entScheme::FVScheme, up::Array{T}) where {T<:AbstractFloat} = OptAlphaGSol{T}(alpha, G, entScheme, up)

# # Constructeur pour OptAlphaGSol
# initialize_OptAlphaGSol()

function divideAlphaG(alphaG::AbstractArray{T}) where {T}
    L = size(alphaG)[1]
    vcat(alphaG[1:Int(L / 2), :], alphaG[1,]), vcat(alphaG[Int(L / 2)+1:end, :], alphaG[Int(L / 2)+1, :])
end

function repasteAlphaG(alpha::AbstractArray{T}, G::AbstractArray{T}) where {T}
    vcat(alpha[begin:end-1], G[begin:end-1])
end

function compute_modified_bounds_alpha(alpha::AbstractArray{T}, sol::OptForEntropySol, alphaScheme::FVScheme) where {T}

    @timeit to "modified bounds computation" begin

        equation, scheme, domain, modifiedDataType, boundsType = sol.equation, sol.scheme, sol.domain, sol.modifiedDataType, sol.boundsType
        Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
        lambd = dt / dx
        sL::Int, sR::Int = get_sL(scheme), get_sR(scheme)
        u, up = sol.u_approx[end-1], sol.u_approx[end]
        # source
        z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx, 1))

        M_delta, m_delta = zeros(T, Nx + 1), zeros(T, Nx + 1)
        for j in 1:Nx
            # Calcul of Bg delta
            ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)

            #uh = scheme_step(equation.source, ut, dt, domain, equation, alphaScheme)
            uh = localized_scheme_step(j, equation.source, ut, dt, domain, equation, alphaScheme)

            zt = isnothing(domain.sourceVec) ? zeros(T, size(ut)) : compute_u_tilde(modifiedDataType, z, j, sL, sR)

            #         ###### This part needs to be completed for equation with source term
            # Dzt = isnothing(domain.DSourceVec) ? zeros(T, size(ut)) : domain.DSourceVec
            # uh = copy(ut)
            # #uh = zeros(T, size(ut))
            # sourceVec = isnothing(sourceTerm(equation, scheme, domain, ut; z=zt, Dz=Dzt)) ? zero(T, Nx) : sourceTerm(equation, scheme, domain, ut; z=zt, Dz=Dzt)
            # #sourceVec = sourceTerm(equation, scheme, domain, ut; z=zt, Dz=Dzt)

            # for k in j-sL-sR+1:j+sR+sL
            #     uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] 
            #     @show typeof(alpha[mod1(k, Nx)])
            #     # @show typeof(uh[mod1(k, Nx), :])
            #     @show typeof(giveNumFlux(alphaScheme.scheme1, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)]))
            #     #uh[mod1(k, Nx), :] .+= alpha[mod1(k, Nx)]
            #     uh[mod1(k, Nx), :] .-= (giveNumFlux(alphaScheme.scheme1, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            #     .-
            #         giveNumFlux(alphaScheme.scheme1, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])).* lambd .* (1-alpha[mod1(k, Nx)])
            #     uh[mod1(k, Nx), :] .-= (
            #             giveNumFlux(alphaScheme.scheme2, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            #             .-
            #             giveNumFlux(alphaScheme.scheme2, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])).* lambd .*(alpha[mod1(k, Nx)])
            #             .+ alpha[mod1(k, Nx)]

            #     uh[mod1(k, Nx), :] .+= dt * sourceVec[mod1(k, Nx),:]
            #         #.-.+ dt * sourceVec[mod1(k, Nx),:] 

            #         # uh[mod1(k, Nx), :] .+=
            #         # -lambd*(1-alpha[mod1(k, Nx)])*(
            #         #     giveNumFlux(alphaScheme.scheme2, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            #         #     .-
            #         #     giveNumFlux(alphaScheme.scheme2, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)]))
            # end

            m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
            m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

            m_delta[j+1], M_delta[j+1] = m[1], M[1]
        end
        m_delta[1], M_delta[1] = m_delta[end], M_delta[end]
        m_delta, M_delta
    end
end

function constraint_alphaG(alphaG::AbstractArray{T}, sol::OptForEntropySol, entScheme::FVScheme) where {T}

    @timeit to "constraint function" begin

        #@code_warntype(divideAlphaG(alphaG))
        alpha, G = divideAlphaG(alphaG)
        equation, domain = sol.equation, sol.domain
        Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
        lambd = dt / dx
        u = sol.u_approx[end-1]
        cons = zero(dx)

        @assert sol.scheme.spaceScheme.CFL_factor == entScheme.spaceScheme.CFL_factor
        alphaScheme = MixedScheme(CFL_factor, sol.scheme, entScheme, alpha)

        # First constraint: Negative diffusion
        #@code_warntype scheme_step(equation.source, u, dt, domain, equation, alphaScheme)
        up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaScheme)
        # bu_delta = zeros(T, Nx)
        # bg_delta = zeros(T, Nx)
        # cons=0.0
        # for j in 1:Nx
        #     bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        #     bg_delta[j] = lambd*((G[j+1] - sol.Gopt[j+1]) - (G[j] - sol.Gopt[j]))
        #     cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
        # end
        #@code_warntype diffusion(u, up_mod, G, dx, dt, equation, domain)
        #@show typeof(G)
        D_mod = diffusion(u, up_mod, G, dx, dt, equation, domain)
        for d in D_mod
            cons += max(0, d)^2
        end
        # display(plot(sol.domain.x, D_mod))
        # plot(sol.domain.x, u)
        # display(plot!(sol.domain.x, up_mod))
        # println("After diffusion, cons="*string(cons))

        # Second constraint: Consistency of numerical entropy flux
        m_delta, M_delta = compute_modified_bounds_alpha(alpha, sol, alphaScheme)
        # @code_warntype compute_modified_bounds_alpha(alpha, sol, alphaScheme)
        #compute_modified_bounds_alpha(alpha, sol, alphaScheme)

        #m_delta, M_delta = compute_G_bounds(u, Nx, dx, dt, equation, domain, alphaScheme, sol.modifiedDataType, sol.boundsType)
        # plot(domain.interfaces, m_delta, label="m")
        # plot!(domain.interfaces, G, label="G")
        # plot!(domain.interfaces, M_delta, label="M")
        # display(title!("bounds of modified scheme"))
        for j in 1:Nx+1
            cons += (dt / dx)^2 * (max(0, m_delta[j] - G[j])^2 + max(0, G[j] - M_delta[j])^2)
        end
        # println("After consistency, cons="*string(cons))


        # Third constraint: alpha must be between specified bounds
        for a in alpha
            cons += max(0, a - 10)^2 + max(0, -9 - a)^2
            #cons += max(0, a-1)^2 + max(0, -a)^2
        end
        # println("After alpha domain, cons="*string(cons))

        # Fourth constraint: maximum principle
        min_u, max_u = minimum(u), maximum(u)
        for j in 1:Nx
            cons += max(0, min_u - up_mod[j]) .^ 2 + max(0, up_mod[j] - max_u) .^ 2
        end


        cons
    end
end

function cost_alphaG(alphaG::AbstractArray{T}, sol::OptForEntropySol, entScheme::FVScheme) where {T}

    @timeit to "cost function" begin

        alpha, G = divideAlphaG(alphaG)
        #sum(alpha[1:end-1].^2)
        alphaScheme = FVScheme(Euler(), MixedScheme(sol.scheme.spaceScheme.CFL_factor, sol.scheme.spaceScheme, entScheme.spaceScheme, alpha))
        u = sol.u_approx[end-1]
        dt = sol.dt_vec[end]
        up_mod = scheme_step(sol.equation.source, u, dt, sol.domain, sol.equation, alphaScheme)

        D = diffusion(u, up_mod, G, sol.domain.dx, dt, sol.equation, sol.domain)

        sum(D .^ 2)
    end

end

# INPUT
# Gent is an array, the numerical entropy flux associated with the entropic scheme

function find_optimal_alphaG(sol::OptForEntropySol, entScheme::FVScheme, Gent)

    @timeit to "find optimal alpha G" begin

        alphaG_init, _p = vcat(zeros(Nx) .+ 1, Gent[1:Nx]), 0.0
        # alphaG_init, _p = vcat(zero(sol.Gopt[:end-1]).+1, Gent[:end-1]), 0.0
        alphaG_init = reshape(alphaG_init, size(alphaG_init)[1])

        cons(res, alphaG, p) = (res .= [constraint_alphaG(alphaG, sol, entScheme)])
        #@show c_init = constraint_alphaG(alphaG_init, sol,  entScheme)
        #@code_warntype constraint_alphaG(alphaG_init, sol,  entScheme)
        cost(alphaG) = cost_alphaG(alphaG, sol, entScheme)

        optprob = OptimizationFunction((x, p) -> cost(x), Optimization.AutoForwardDiff(); cons=cons)
        prob = @timeit to "defining" OptimizationProblem(optprob, alphaG_init, _p, lcons=[-1000.0], ucons=[0.0])
        solAlphaG = @timeit to "solving" solve(prob, Ipopt.Optimizer(); max_iter=100, tol=1e-2)

        alphaSol, GSol = divideAlphaG(solAlphaG.u)
        equation, domain = sol.equation, sol.domain
        dt = sol.dt_vec[end]
        u = sol.u_approx[end-1]
        alphaScheme = MixedScheme(CFL_factor, sol.scheme, entScheme, alphaSol)
        up_mod = scheme_step(equation.source, u, dt, domain, equation, alphaScheme)

        OptAlphaGSol(alphaSol, GSol, entScheme, up_mod)
    end
end