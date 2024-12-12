abstract type OptimFunctional end
struct SquareMinFun <: OptimFunctional end
struct AbsMinFun <: OptimFunctional end
struct SqrtMinFun <: OptimFunctional end
struct HeavyConsistencyMinFun <:OptimFunctional end
struct WeightedMinFun <: OptimFunctional end

abstract type InitGuess end
struct MeanInitGuess <: InitGuess end
struct NullInitGuess <: InitGuess end

function compute_u_tilde(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int)
    
    K = computeK(KFun, extractLocalData(u, j, sL, sR))
    Nx= size(u)[1]
    ut = zero(u)
    
    for i in 1:Nx
        if (j - sL + 1 <= i && i <= j + sR) || (((j - sL + 1) % Nx <= i) && (j - sL + 1 < 1)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
            ut[i,:] = u[i,:]
        else
            ut[i,:] .= K
        end
    end

    ut

end


function compute_u_tilde(::AsymmetricModifiedData, u, j::Int, sL::Int, sR::Int)

    Nx = size(u)[1]
    ut = zero(u)

    ileft, iright = j - sL + 1, j + sR
    jleft, jright = j - sL - sR + 2, j + sL + sR - 1

    for k in 1:Nx

        if (
            (ileft <= k && iright >= k)
            || (ileft % Nx <= k && ileft < 1)
            || (iright % Nx >= k && iright > Nx)
        )
            ut[k,:] = u[k,:]
        elseif (jleft <= k && ileft > k) || (
            jleft % Nx <= k && jleft < 1 && k < ileft % Nx
        )
            ut[k,:] = u[mod1(ileft, Nx),:]
        elseif (iright < k && jright >= k) || (
            jright % Nx >= k && jright > Nx && k > iright % Nx
        )
            ut[k,:] = u[mod1(iright, Nx),:]
        elseif jright < k
            ut[k,:] = u[mod1(iright, Nx),:]
        elseif jleft > k
            ut[k,:] = u[mod1(ileft, Nx),:]
        end
    end
    ut
end

compute_z_tilde(::NullSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR) = nothing
function compute_z_tilde(zbSource::ZbSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR)
    z = reshape(domain.sourceVec, (domain.Nx, 1))
    compute_u_tilde(modifiedDataType, z, j, sL, sR)
end

function compute_u_hat(ns::NullSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=nothing)

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    # numericalFluxMat = giveNumFlux(zb, method, equation, v; domain=domain)

    # for k in j-sL-sR+1:j+sR+sL

    #     uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
    #         numFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :])
    #         .-
    #         numFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :])
    #     )

    # end
    for k in j-sL-sR+1:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            numFlux(method, equation, extract_data_stencil(equation, ut, mod1(k, Nx), sL, sR)...)
            .-
            numFlux(method, equation, extract_data_stencil(equation, ut, mod1(k-1, Nx), sL, sR)...)
        )

    end
    #scheme_step(ns, ut, dt, domain, equation, method)
    return uh
end

function compute_u_hat(::ZbSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=zero(ut))

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    sourceVec = sourceTerm(equation, method, domain, ut; z=zt)

    for k in j-sL-sR+1:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            giveNumFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            .-
            giveNumFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])) .+ dt * sourceVec[mod1(k, Nx),:]
        

    end

    return uh
end

function initBounds(KFun::SymmetricModifiedData, equation::Equation, u, j, sL, sR, z=nothing)
    GK = get_G(equation, computeK(KFun, extractLocalData(u, j, sL, sR)), computeZ(KFun, z, j, sL, sR))
    return GK[1], GK[1]
end

function initBounds(::AsymmetricModifiedData, equation::Equation, u, j, sL, sR, z=nothing)
    Nx = size(u)[1]
    zL = isnothing(z) ? nothing : z[mod1(j + sR, Nx),:]
    zR = isnothing(z) ? nothing : z[mod1(j - sL + 1, Nx),:]
    Gm = get_G(equation, u[mod1(j + sR, Nx),:], zR)
    GM = get_G(equation, u[mod1(j - sL + 1, length(u)),:], zL)
    Gm[1], GM[1]
end

function updateBounds!(::SymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    for k in j-sL-sR+1:j
        M = M .+ dx / dt .* (get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]))
    end

    for k in j+1:j+sL+sR
        m = m .+ dx / dt .* (get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]))
    end
    
    m, M

end

function updateBounds!(::AsymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    for k in j-sL-sR+2:j
        M = M .+ dx / dt .* (get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]))
    end

    for k in j+1:j+sL+sR-1
        m = m .+ dx / dt .* (get_eta(equation, uh[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]))
    end
    
    m, M

end

function updateBounds!(KFun::SymmetricModifiedData, ::LightBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    Deta_uj = get_D_eta(equation, u[j,:], zt[j,:])
    eta_uj = get_eta(equation, u[j,:], zt[j,:])
    fj = giveNumFlux(method, equation, u[j,:], u[mod1(j+1, Nx),:]; zL=zt[j], zR=zt[mod1(j + 1, Nx)])
    K, Z = computeK(KFun, extractLocalData(u, j, sL, sR)), computeZ(KFun, z, j, sL, sR)

    M = M .+ Deta_uj * (fj .- get_flux(equation, K, Z))
    M = M .+ dx/dt * (j + 1 - (j - sL - sR + 1))* (Deta_uj.*u[j,:] .- eta_uj)

    m = m .+ Deta_uj * (fj .- get_flux(equation, K, Z))
    m = m .+ dx/dt * ((j + sL + sR + 1) - j - 1)* (eta_uj .- Deta_uj.*u[j,:])
    
    for k in j-sL-sR+1:j
        M = M .+ dx / dt .* (get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]) .- Deta_uj * ut[mod1(k, Nx),:])
    end

    for k in j+1:j+sL+sR
        m = m .+ dx / dt .* (Deta_uj * ut[mod1(k, Nx),:] .- get_eta(equation, ut[mod1(k, Nx),:], zt[mod1(k, Nx),:]))
        
    end
    
    m, M

end

function compute_G_bounds(u, Nx, dx, dt, equation::Equation, domain::Domain, method::FVMethod, modifiedDataType::ModifiedDataType, boundsType::BoundsType=NormalBounds())

    M_vec, m_vec = zeros(Nx + 1), zeros(Nx + 1)
    sL, sR = get_sL(method), get_sR(method)

    # source
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    for j in 1:Nx
        ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
        uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)
        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        m_vec[j+1], M_vec[j+1] = m[1], M[1]

    end

    m_vec[1], M_vec[1] = m_vec[end], M_vec[end]

    m_vec, M_vec

    
end

function J(::SquareMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)
    JD = zero(eltype(u))
    JC = zero(eltype(u))
    z = isnothing(domain.sourceVec) ? zeros(Nx,1) : domain.sourceVec

    for j in 1:Nx
        JD += max(0, get_eta(equation, up[j,:], z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2
    end

    for j in 1:Nx+1
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end

    JD + JC + 1
end

function J(::HeavyConsistencyMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = domain.sourceVec

    for j in 1:Nx
        JD += max(0, get_eta(equation, up[j,:], z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2
    end

    for j in 1:Nx+1
        JC += (dt / dx) * (max(0, (gamma[j] - M_vec[j])) + max(0, 1 * (m_vec[j] - gamma[j])))
    end

    JD + JC + 1
end

function J(::SqrtMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = domain.sourceVec

    for j in 1:Nx
        JD += sqrt(max(0, get_eta(equation, up[j,:], z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j])))
    end

    for j in 1:Nx+1
        JC += sqrt((dt / dx)) * (sqrt(max(0, (gamma[j] - M_vec[j]))) + sqrt(max(0, 1 * (m_vec[j] - gamma[j]))))
    end

    JD + JC
end


function J(::AbsMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = domain.sourceVec

    for j in 1:Nx
        JD += max(0, get_eta(equation, up[j,:], z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))
    end

    for j in 1:Nx+1
        JC += (dt / dx) * (max(0, (gamma[j] - M_vec[j])) + max(0, 1 * (m_vec[j] - gamma[j])))
    end

    JD + JC
end

function J(::WeightedMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = isnothing(domain.sourceVec) ? zeros(Nx,1) : domain.sourceVec

    for j in 1:Nx
        JD += max(0, get_eta(equation, up[j,:], z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2 * ((M_vec[j+1] - m_vec[j+1]) >=0)
        #JD += max(0, get_eta(equation, up[j,:]; z=z[j])[1] - get_eta(equation, u[j,:]; z=z[j])[1] + dt / dx * (gamma[j+1] - gamma[j])) * ((M_vec[j+1] - m_vec[j+1]) >=0)
    end

    for j in 1:Nx+1
        #@show (M_vec[j] - m_vec[j]) >=0
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2) * ((M_vec[j] - m_vec[j]) >=0)
        #JC += (dt / dx) * (max(0, (gamma[j] - M_vec[j])) + max(0, 1 * (m_vec[j] - gamma[j]))) * ((M_vec[j] - m_vec[j]) >=0)
    end

    JD + JC + 1
end

function diffusion(u, up, gamma::AbstractArray{U}, dx::T, dt::T, equation::Equation, domain::Domain{T}) where T<:Real where U
    # z = isnothing(domain.sourceVec) ? zeros(T, (Nx,1)) : domain.sourceVec
    z = manageSource(domain)
    #@code_warntype(manageSource(domain))
    [get_eta(equation, up[i,:], z[i])[1] - get_eta(equation, u[i,:], z[i])[1] for i in 1:length(u[:,1])] .+ dt / dx * (gamma[2:end] - gamma[1:end-1])
end

function consistency(::SquareMinFun, gamma, Nx, dx, dt, m_vec, M_vec)
    consistency_vec = zeros(Nx+1)
    for j in 1:Nx+1
        consistency_vec[j] = (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end
    consistency_vec
end
function consistency(::WeightedMinFun, gamma, Nx, dx, dt, m_vec, M_vec)
    consistency_vec = zeros(Nx+1)
    for j in 1:Nx+1
        consistency_vec[j] = (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2) * ((M_vec[j] - m_vec[j]) >=0)
    end
    consistency_vec
end
consistency(::OptimFunctional, gamma, Nx, dx, dt, m_vec, M_vec) = @warn "Missing consistency function"

initial_guess(::MeanInitGuess, m_vec, M_vec) = 0.5 * (m_vec + M_vec)
initial_guess(::NullInitGuess, m_vec, M_vec) = zero(m_vec)

function optimize_for_entropy(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(get_sL(method), get_sR(method)), boundsType::BoundsType=NormalBounds(), optimFunctional::OptimFunctional=SquareMinFun(), initGuess::InitGuess=MeanInitGuess(), kwargs...)

    Nx, dx = domain.Nx, domain.dx
    FVsol = fv_solve(domain, u_init, equation, method)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, method, modifiedDataType, boundsType)
    gamma_init = initial_guess(initGuess, m_vec, M_vec)

    u, up, dt = u_approx[end-1], u_approx[end], dt_vec[end]
    sol = optimize(gamma -> J(optimFunctional, gamma, u_approx[end-1], u_approx[end], Nx, dx, dt_vec[end], m_vec, M_vec, equation, domain), gamma_init; kwargs...)#g_tol=1e-10, iterations=100000)#; g_tol=1e-10)#; autodiff=:forward)#, kwargs...)
    Gopt, Jopt = Optim.minimizer(sol), Optim.minimum(sol)

    Dopt = diffusion(u_approx[end-1], u_approx[end], Gopt, dx, dt_vec[end], equation, domain)
    Copt = consistency(optimFunctional, Gopt, Nx, dx, dt_vec[end], m_vec, M_vec)
    # @show typeof(u_approx)
    # @show typeof(domain.dx)
    OptForEntropySol(domain, equation, method, modifiedDataType, boundsType, u_approx, dt_vec, Gopt, Jopt, Dopt, Copt, m_vec, M_vec, sol, "")

end