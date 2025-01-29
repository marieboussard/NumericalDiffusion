abstract type OptimFunctional end
struct SquareMinFun <: OptimFunctional end
struct AbsMinFun <: OptimFunctional end
struct SqrtMinFun <: OptimFunctional end
struct HeavyConsistencyMinFun <: OptimFunctional end
struct WeightedMinFun <: OptimFunctional end

abstract type InitGuess end
struct MeanInitGuess <: InitGuess end
struct NullInitGuess <: InitGuess end

function compute_u_tilde(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int)

    K = computeK(KFun, extractLocalData(u, j, sL, sR))
    Nx = size(u)[1]
    ut = zero(u)

    for i in 1:Nx
        if (j - sL + 1 <= i && i <= j + sR) || (((j - sL + 1) % Nx <= i) && (j - sL + 1 < 1)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
            ut[i, :] = u[i, :]
        else
            ut[i, :] .= K
        end
    end

    ut

end

function compute_u_tilde(::AsymmetricModifiedData, u, j::Int, sL::Int, sR::Int)
    Nx = size(u)[1]
    u_shift, ut_shift, ut = zero(u), zero(u), zero(u)

    #@show j-sL-sR+2
    # First the data is shifted so that the index j-sL-sR+2 is brought to position 1
    for k in 1:Nx
        u_shift[k] = u[mod1(k + j - 2 * sL - sR + 1, Nx)]
    end
    #@show u_shift

    # Then, the shifted data is modified
    for k in 1:sL+sR-1
        ut_shift[k] = u_shift[sL+sR]
    end
    for k in sL+sR:2*sR+2*sL-1
        ut_shift[k] = u_shift[k]
    end
    for k in 2*sR+2*sL:3*sL+3*sR-2
        ut_shift[k] = u_shift[2*sR+2*sL-1]
    end

    #@show ut_shift

    # Finally, we apply the inverse shifting
    for k in 1:Nx
        ut[k] = ut_shift[mod1(k - j + 2 * sL + sR - 1, Nx)]
    end

    #@show ut

    ut

end

# function compute_u_tilde(::AsymmetricModifiedData, u, j::Int, sL::Int, sR::Int)

#     Nx = size(u)[1]
#     ut = zero(u)

#     ileft, iright = j - sL + 1, j + sR
#     jleft, jright = j - sL - sR + 2, j + sL + sR - 1

#     for k in 1:Nx

#         if (
#             (ileft <= k && iright >= k)
#             || (ileft % Nx <= k && ileft < 1)
#             || (iright % Nx >= k && iright > Nx)
#         )
#             ut[k,:] = u[k,:]
#         elseif (jleft <= k && ileft > k) || (
#             jleft % Nx <= k && jleft < 1 && k < ileft % Nx
#         )
#             ut[k,:] = u[mod1(ileft, Nx),:]
#         elseif (iright < k && jright >= k) || (
#             jright % Nx >= k && jright > Nx && k > iright % Nx
#         )
#             ut[k,:] = u[mod1(iright, Nx),:]
#         elseif jright < k
#             ut[k,:] = u[mod1(iright, Nx),:]
#         elseif jleft > k
#             ut[k,:] = u[mod1(ileft, Nx),:]
#         end
#     end
#     ut
# end

compute_z_tilde(::NullSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR) = nothing
function compute_z_tilde(zbSource::ZbSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR)
    z = reshape(domain.sourceVec, (domain.Nx, 1))
    compute_u_tilde(modifiedDataType, z, j, sL, sR)
end

function compute_u_hat(ns::NullSource, ut, dx, dt, j, domain::Domain, equation::Equation, scheme::FVScheme; zt=nothing)

    uh = copy(ut)
    Nx, p = size(ut)
    sL, sR = get_sL(scheme), get_sR(scheme)

    # numericalFluxMat = giveNumFlux(zb, method, equation, v; domain=domain)

    # for k in j-sL-sR+1:j+sR+sL

    #     uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
    #         numFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :])
    #         .-
    #         numFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :])
    #     )

    # end

    ######### TO DO: CHANGE TO IMPROVE COMPUTATIONNAL TIME ##################################
    # for k in j-sL-sR+1:j+sR+sL

    #     uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
    #         numFlux(scheme, equation, extract_data_stencil(equation, ut, mod1(k, Nx), sL, sR)...)
    #         .-
    #         numFlux(scheme, equation, extract_data_stencil(equation, ut, mod1(k-1, Nx), sL, sR)...)
    #     )

    # end
    # return uh
    scheme_step(ns, ut, dt, domain, equation, scheme)

end

function compute_u_hat(zb::ZbSource, ut, dx, dt, j, domain::Domain, equation::Equation, scheme::FVScheme; zt=zero(ut))

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(scheme), get_sR(scheme)

    sourceVec = sourceTerm(equation, scheme.spaceScheme, domain, ut; z=zt)

    # for k in j-sL-sR+1:j+sR+sL

    # uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
    #     giveNumFlux(scheme, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
    #     .-
    #     giveNumFlux(scheme, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])) .+ dt * sourceVec[mod1(k, Nx),:]


    # end

    # return uh
    scheme_step(zb, ut, dt, domain, equation, scheme; z=zt)
end

function initBounds(KFun::SymmetricModifiedData, equation::Equation, u, j, sL, sR, z=nothing)
    GK = G(equation, computeK(KFun, extractLocalData(u, j, sL, sR)), computeZ(KFun, z, j, sL, sR))
    return GK[1], GK[1]
end

function initBounds(::AsymmetricModifiedData, equation::Equation, u, j, sL, sR, z=nothing)
    Nx = size(u)[1]
    zL = isnothing(z) ? nothing : z[mod1(j + sR, Nx), :]
    zR = isnothing(z) ? nothing : z[mod1(j - sL + 1, Nx), :]
    Gm = G(equation, u[mod1(j + sR, Nx), :], zR)
    GM = G(equation, u[mod1(j - sL + 1, length(u)), :], zL)
    Gm[1], GM[1]
end

function updateBounds!(::SymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    for k in j-sL-sR+1:j
        M = M .+ dx / dt .* (eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]) .- eta(equation, uh[mod1(k, Nx), :], zt[mod1(k, Nx), :]))
    end

    for k in j+1:j+sL+sR
        m = m .+ dx / dt .* (eta(equation, uh[mod1(k, Nx), :], zt[mod1(k, Nx), :]) .- eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]))
    end

    m, M

end

function updateBounds!(::AsymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    for k in j-sL-sR+2:j
        M = M .+ dx / dt .* (eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]) .- eta(equation, uh[mod1(k, Nx), :], zt[mod1(k, Nx), :]))
    end

    for k in j+1:j+sL+sR-1
        m = m .+ dx / dt .* (eta(equation, uh[mod1(k, Nx), :], zt[mod1(k, Nx), :]) .- eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]))
    end

    m, M

end

function updateBounds!(KFun::SymmetricModifiedData, ::LightBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    Deta_uj = D_eta(equation, ut[j, :], zt[j, :])
    eta_uj = eta(equation, ut[j, :], zt[j, :])
    #fj = giveNumFlux(scheme, equation, ut[j, :], ut[mod1(j + 1, Nx), :]; zL=zt[j], zR=zt[mod1(j + 1, Nx)])
    fj = numFlux(scheme, equation, extract_data_stencil(ut, j, sL, sR); z = extract_data_stencil(zt, j, sL, sR), dt=dt, domain=domain)
    K, Z = computeK(KFun, extractLocalData(ut, j, sL, sR)), computeZ(KFun, zt, j, sL, sR)

    M = M .+ Deta_uj .* (fj .- flux(equation, K, Z))
    M = M .+ dx / dt .* (j + 1 - (j - sL - sR + 1)) * (Deta_uj .* ut[j, :] .- eta_uj)

    m = m .+ Deta_uj .* (fj .- flux(equation, K, Z))
    m = m .+ dx / dt .* ((j + sL + sR + 1) - j - 1) * (eta_uj .- Deta_uj .* ut[j, :])

    for k in j-sL-sR+1:j
        M = M .+ dx / dt .* (eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]) .- Deta_uj .* ut[mod1(k, Nx), :])
    end

    for k in j+1:j+sL+sR
        m = m .+ dx / dt .* (Deta_uj .* ut[mod1(k, Nx), :] .- eta(equation, ut[mod1(k, Nx), :], zt[mod1(k, Nx), :]))

    end

    m, M

end

function compute_local_bounds(boundsType::BoundsType, u, dx, dt, j, sL, sR, z, domain::Domain, equation::Equation, scheme::FVScheme, modifiedDataType::ModifiedDataType)
    ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
    zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
    uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, scheme; zt=zt)
    m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
    m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)
end

function compute_local_bounds(::SimpleBounds, u, dx, dt, j, sL, sR, z, domain::Domain, equation::Equation, scheme::FVScheme, args...)
    F_j = numFlux(scheme, equation, extract_data_stencil(u, j, sL, sR); dt=dt, domain=domain)
    uh_j = u[j, :] .- dt / dx * (F_j .- flux(equation, u[j, :])[1])
    uh_jp = u[mod1(j + 1, Nx), :] .- dt / dx * (flux(equation, u[mod1(j + 1, Nx), :])[1] .- F_j)
    m = G(equation, u[mod1(j + 1, Nx), :], z[mod1(j + 1, Nx)])[1] + dx / dt * (eta(equation, uh_jp, z[mod1(j + 1, Nx)])[1] - eta(equation, u[mod1(j + 1, Nx), :], z[mod1(j + 1, Nx)])[1])
    M = G(equation, u[j, :], 0.0)[1] - dx / dt * (eta(equation, uh_j, 0.0)[1] - eta(equation, u[j, :], 0.0)[1])
    m, M
end

# Removing ambiguity
compute_local_bounds(::SimpleBounds, u, dx, dt, j, sL, sR, z, domain::Domain, equation::Equation, scheme::FVScheme, modifiedDataType::ModifiedDataType) = compute_local_bounds(SimpleBounds(), u, dx, dt, j, sL, sR, z, domain, equation, scheme)

function compute_G_bounds(u, Nx, dx, dt, equation::Equation, domain::Domain, scheme::FVScheme, modifiedDataType::ModifiedDataType, boundsType::BoundsType=NormalBounds())

    M_vec, m_vec = zeros(Nx + 1), zeros(Nx + 1)
    sL, sR = get_sL(scheme), get_sR(scheme)

    # source
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx, 1))

    for j in 1:Nx
        # ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
        # uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, scheme; zt=zt)
        # m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)
        # m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)
        m, M = compute_local_bounds(boundsType, u, dx, dt, j, sL, sR, z, domain, equation, scheme, modifiedDataType)

        m_vec[j+1], M_vec[j+1] = m[1], M[1]

    end

    m_vec[1], M_vec[1] = m_vec[end], M_vec[end]

    m_vec, M_vec


end

function J(::SquareMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)
    JD = zero(eltype(u))
    JC = zero(eltype(u))
    z = isnothing(domain.sourceVec) ? zeros(Nx, 1) : domain.sourceVec

    for j in 1:Nx
        JD += max(0, eta(equation, up[j, :], z[j])[1] - eta(equation, u[j, :], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2
    end

    for j in 1:Nx+1
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end

    JD + JC# + 1
end

function J(::HeavyConsistencyMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = domain.sourceVec

    for j in 1:Nx
        JD += max(0, eta(equation, up[j, :], z[j])[1] - eta(equation, u[j, :], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2
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
        JD += sqrt(max(0, eta(equation, up[j, :], z[j])[1] - eta(equation, u[j, :], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j])))
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
        JD += max(0, eta(equation, up[j, :], z[j])[1] - eta(equation, u[j, :], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))
    end

    for j in 1:Nx+1
        JC += (dt / dx) * (max(0, (gamma[j] - M_vec[j])) + max(0, 1 * (m_vec[j] - gamma[j])))
    end

    JD + JC
end

function J(::WeightedMinFun, gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = isnothing(domain.sourceVec) ? zeros(Nx, 1) : domain.sourceVec

    for j in 1:Nx
        JD += max(0, eta(equation, up[j, :], z[j])[1] - eta(equation, u[j, :], z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2 * ((M_vec[j+1] - m_vec[j+1]) >= 0)
        #JD += max(0, eta(equation, up[j,:]; z=z[j])[1] - eta(equation, u[j,:]; z=z[j])[1] + dt / dx * (gamma[j+1] - gamma[j])) * ((M_vec[j+1] - m_vec[j+1]) >=0)
    end

    for j in 1:Nx+1
        #@show (M_vec[j] - m_vec[j]) >=0
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2) * ((M_vec[j] - m_vec[j]) >= 0)
        #JC += (dt / dx) * (max(0, (gamma[j] - M_vec[j])) + max(0, 1 * (m_vec[j] - gamma[j]))) * ((M_vec[j] - m_vec[j]) >=0)
    end

    JD + JC + 1
end

function diffusion(u, up, gamma::AbstractArray{U}, dx::T, dt::T, equation::Equation, domain::Domain{T}) where {T<:Real} where {U}
    # z = isnothing(domain.sourceVec) ? zeros(T, (Nx,1)) : domain.sourceVec
    z = manageSource(domain)
    #@code_warntype(manageSource(domain))
    [eta(equation, up[i, :], z[i])[1] - eta(equation, u[i, :], z[i])[1] for i in 1:length(u[:, 1])] .+ dt / dx * (gamma[2:end] - gamma[1:end-1])
end

diffusion(sol::OptForEntropySol, Gexact) = diffusion(sol.u_approx[end-1], sol.u_approx[end], Gexact, sol.domain.dx, sol.dt_vec[end], sol.equation, sol.domain)

function consistency(::SquareMinFun, gamma, Nx, dx, dt, m_vec, M_vec)
    consistency_vec = zeros(Nx + 1)
    for j in 1:Nx+1
        consistency_vec[j] = (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end
    consistency_vec
end
function consistency(::WeightedMinFun, gamma, Nx, dx, dt, m_vec, M_vec)
    consistency_vec = zeros(Nx + 1)
    for j in 1:Nx+1
        consistency_vec[j] = (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2) * ((M_vec[j] - m_vec[j]) >= 0)
    end
    consistency_vec
end
consistency(::OptimFunctional, gamma, Nx, dx, dt, m_vec, M_vec) = @warn "Missing consistency function"

initial_guess(::MeanInitGuess, m_vec, M_vec) = 0.5 * (m_vec + M_vec)
initial_guess(::NullInitGuess, m_vec, M_vec) = zero(m_vec)

function optimize_for_entropy(u_init, domain::Domain, equation::Equation, scheme::FVScheme; modifiedDataType::ModifiedDataType=meanK(get_sL(scheme), get_sR(scheme)), boundsType::BoundsType=NormalBounds(), optimFunctional::OptimFunctional=SquareMinFun(), initGuess::InitGuess=MeanInitGuess(), vectorInitGuess=nothing, kwargs...)

    Nx, dx = domain.Nx, domain.dx
    FVsol = fv_solve(domain, u_init, equation, scheme)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, scheme, modifiedDataType, boundsType)
    if !isnothing(vectorInitGuess)
        gamma_init = vectorInitGuess
    else
        gamma_init = initial_guess(initGuess, m_vec, M_vec)
    end

    u, up, dt = u_approx[end-1], u_approx[end], dt_vec[end]
    sol = optimize(gamma -> J(optimFunctional, gamma, u_approx[end-1], u_approx[end], Nx, dx, dt_vec[end], m_vec, M_vec, equation, domain), gamma_init; g_tol=1e-20, iterations=200000, method=LBFGS(), autodiff=:forward, kwargs...)#g_tol=1e-10, iterations=100000)#; g_tol=1e-10)#; autodiff=:forward)#, kwargs...),

    Gopt, Jopt = Optim.minimizer(sol), Optim.minimum(sol)

    Dopt = diffusion(u_approx[end-1], u_approx[end], Gopt, dx, dt_vec[end], equation, domain)
    Copt = consistency(optimFunctional, Gopt, Nx, dx, dt_vec[end], m_vec, M_vec)

    OptForEntropySol(sol, domain, equation, scheme, modifiedDataType, boundsType, u_approx, dt_vec, Gopt, Jopt, Dopt, Copt, m_vec, M_vec, sol, "")

end