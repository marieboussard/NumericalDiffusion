abstract type GmidShape end
struct ConsistencyMeanShape <: GmidShape end
struct RusanovLikeShape <: GmidShape end

function compute_modified_bounds_f_delta(f_delta::AbstractArray{T}, sol::OptForEntropySol) where T
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
                .-
                lambd*(f_delta[mod1(k, Nx)] - f_delta[mod1(k-1, Nx)])
        end

        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        m_delta[j+1], M_delta[j+1] = m[1], M[1]
    end
    m_delta, M_delta
end

function compute_G_mid_f_delta(f_delta::AbstractArray{T}, sol::OptForEntropySol, ::ConsistencyMeanShape) where T
    m_delta, M_delta = compute_modified_bounds_f_delta(f_delta, sol)
    0.5*(m_delta.+M_delta)
end

function compute_G_mid(A::AbstractArray{T}, sol::OptForEntropySol, ::ConsistencyMeanShape) where T
    compute_G_mid_f_delta(f_delta_from_A(A, sol), sol, ConsistencyMeanShape())
end

function modified_scheme_step_f_delta(f_delta::AbstractArray{T}, sol::OptForEntropySol) where T
    equation, method, domain = sol.equation, sol.method, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    u = sol.u_approx[end-1]
    # source
    z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))
    
    p = get_unknowns_number(equation)
    numericalFluxMat = zeros(T, (Nx+1, p))
    up = zeros(T, (Nx, p))

    for i ∈ 2:Nx
        numericalFluxMat[i,:] = giveNumFlux(method, equation, u[i-1,:], u[i,:]; zL=z[i-1], zR=z[i])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, u[end,:], u[1,:]; zL=z[end], zR=z[1])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]

    # for j in 1:Nx
    #     up[j] = u[j] - lambd*(numericalFluxMat[j] - numericalFluxMat[mod1(j-1,Nx)]) - lambd*(f_delta[j+1] - f_delta[j])
    # end

    up = u - lambd * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:]) - lambd*(f_delta[2:end,:] - f_delta[1:end-1,:]) #+ dt * sourceTerm(equation, method, domain, v)

    #@show up.-sol.u_approx[end]
    up

end

function modified_scheme_step(A::AbstractArray{T}, sol::OptForEntropySol) where T
    modified_scheme_step_f_delta(f_delta_from_A(A, sol), sol)
end

function compute_G_mid(A::AbstractArray{T}, sol::OptForEntropySol, ::RusanovLikeShape) where T
    
    equation, domain = sol.equation, sol.domain
    Nx = domain.Nx
    u = sol.u_approx[end-1]
    # source
    z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    Gmid = zeros(T, Nx+1)
    for j in 1:Nx-1
        Gmid[j+1] = sol.Gopt[j+1] + 0.5*A[j+1]*(get_eta(equation, u[j+1]; z=z[j+1]) - get_eta(equation, u[j]; z=z[j]))
    end
    Gmid[end] = sol.Gopt[end] + 0.5*A[end]*(get_eta(equation, u[1]; z=z[1]) - get_eta(equation, u[end]; z=z[end]))
    Gmid[1] = Gmid[end]

    #@show Gmid.-sol.Gopt
    Gmid
end

#=
function modified_scheme_step(A::AbstractArray{T}, sol::OptForEntropySol, ::RusanovLikeShape) where T
end
=#

function diffusion_from_A(A::AbstractArray{T}, sol::OptForEntropySol; gmidShape::GmidShape=ConsistencyMeanShape()) where T

    # Extracting data from solution of entropy optimizing process
    domain, equation, u = sol.domain, sol.equation, sol.u_approx[end-1]
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    z = isnothing(domain.sourceVec) ? zeros(T, (Nx,1)) : domain.sourceVec

    # Constructing approximated solution and numerical entropy flux for modified scheme associated with A
    # f_delta = f_delta_from_A(A, sol)
    # G_mid = compute_G_mid_f_delta(f_delta, sol)
    # up = modified_scheme_step_f_delta(f_delta, sol)
    up, Gmid = compute_u_and_Gmid_from_A(A, sol, gmidShape)
    
    # Calculating numerical diffusion
    Dmid = [get_eta(equation, up[i,:]; z=z[i])[1] - get_eta(equation, u[i,:]; z=z[i])[1] for i in 1:length(u[:,1])] + dt / dx * (Gmid[2:end] - Gmid[1:end-1])
    #@show Dmid.-sol.Dopt
end

#=
function constraint_f_delta(f_delta::AbstractArray{T}, sol::OptForEntropySol) where T
    
    equation, method, domain, modifiedDataType, boundsType = sol.equation, sol.method, sol.domain, sol.modifiedDataType, sol.boundsType
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    #sL, sR = get_sL(method), get_sR(method)
    u, up = sol.u_approx[end-1], sol.u_approx[end]
    # source
    #z = isnothing(domain.sourceVec) ? zeros(T, (Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))
    
    #p = get_unknowns_number(equation)
    #numericalFluxMat = zeros(T, (Nx+1, p))

    # Computing solution and numerical entropy flux for the modified scheme
    up_mod = modified_scheme_step_f_delta(f_delta, sol)
    Gmid = compute_G_mid_f_delta(f_delta, sol)
    
    #=
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

        =#

    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    cons=0.0
    for j in 1:Nx
        # Calcul of Bu delta
        # bu_delta[j] = get_eta(equation, u[j] - lambd*(numericalFluxMat[j] - numericalFluxMat[mod1(j-1,Nx)]) - lambd*(f_delta[j+1] - f_delta[j])) - get_eta(equation, up[j])
        bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])

        # Calcul of Bg delta
        # bg_delta[j] = lambd*(((M_delta[j+1]+m_delta[j+1])/2 - sol.Gopt[j+1]) - ((M_delta[j]+m_delta[j])/2 - sol.Gopt[j]))
        bg_delta[j] = lambd*((Gmid[j+1] - sol.Gopt[j+1]) - (Gmid[j] - sol.Gopt[j]))

        cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
        
    end
    cons
    
    #f_delta[1]
end
=#

function f_delta_from_A(A::AbstractArray{T}, sol::OptForEntropySol) where T

    # A should be a vector

    # First case : we use only initial data
    u = sol.u_approx[end-1]

    # # Second case : we use data estimated by non entropic scheme
    # u = sol.u_approx[end]

    Nx, p = size(u)
    f_delta = zeros(T, (Nx+1, p))
    for i in 1:Nx+1
        f_delta[i,:] = -0.5*A[i].* (u[mod1(i+1, Nx), :] - u[mod1(i, Nx), :])
    end
    #zero(A) .+ (1-A[1]^2)
    #@show f_delta
    f_delta
end

#=
function constraint_A(A, sol::OptForEntropySol)
    constraint_f_delta(f_delta_from_A(A, sol), sol)
end
=#

function compute_u_and_Gmid_from_A(A::AbstractArray{T}, sol::OptForEntropySol, gmidShape::GmidShape) where T
    modified_scheme_step(A, sol), compute_G_mid(A, sol, gmidShape)
end

function compute_u_and_Gmid_from_A(A::AbstractArray{T}, sol::OptForEntropySol, ::ConsistencyMeanShape) where T
    f_delta = f_delta_from_A(A, sol)
    modified_scheme_step_f_delta(f_delta, sol), compute_G_mid_f_delta(f_delta, sol, ConsistencyMeanShape())
end

function constraint_A(A::AbstractArray{T}, sol::OptForEntropySol; gmidShape::GmidShape=ConsistencyMeanShape) where T
    equation, domain = sol.equation, sol.domain
    Nx, dx, dt = domain.Nx, domain.dx, sol.dt_vec[end]
    lambd = dt/dx
    up = sol.u_approx[end]

    # Computing solution and numerical entropy flux for the modified scheme
    up_mod, Gmid = compute_u_and_Gmid_from_A(A, sol, gmidShape)
    bu_delta = zeros(T, Nx)
    bg_delta = zeros(T, Nx)
    cons=0.0
    for j in 1:Nx
        # Calcul of Bu delta
        bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        # Calcul of Bg delta
        bg_delta[j] = lambd*((Gmid[j+1] - sol.Gopt[j+1]) - (Gmid[j] - sol.Gopt[j]))
        # Adding the local constraint
        cons += min(0, min(0,-sol.Dopt[j])-bu_delta[j]-bg_delta[j])^2
        
    end
    cons
end

function find_optimal_A(sol::OptForEntropySol; gmidShape::GmidShape=ConsistencyMeanShape())
    @show A_init, _p = zero(sol.Gopt).-100.0, 0.0
    cons(res, A, p) = (res .= [constraint_A(A, sol; gmidShape=gmidShape)])
    @show constraint_A(A_init, sol; gmidShape=gmidShape)

    #optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoForwardDiff(); cons = cons)
    optprob = OptimizationFunction((x,p) -> sum(x.^2), Optimization.AutoForwardDiff(); cons = cons)
    #optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoZygote(); cons = cons)
    #optprob = OptimizationFunction((x,p) -> abs.(x); cons = cons)
    # optprob = OptimizationFunction((x,p) -> abs.(x), Optimization.AutoModelingToolkit(); cons = cons)
    prob = OptimizationProblem(optprob, A_init, _p, lcons = [-1000.0], ucons = [0.0]) # -Inf
    #solA = solve(prob, IPNewton())
    #solA = solve(prob, BFGS())
    solA = solve(prob, Ipopt.Optimizer());
end