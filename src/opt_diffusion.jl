abstract type ModifiedDataType end
abstract type SymmetricModifiedData <: ModifiedDataType end
struct CLModifiedData <: SymmetricModifiedData
    weights
end
struct MaxModifiedData <: SymmetricModifiedData end
struct MinModifiedData <: SymmetricModifiedData end
struct AsymmetricModifiedData <: ModifiedDataType end

function extractLocalData(u, j, sL, sR)

    # To keep only cells that K takes as arguments
    #Nx = length(u)
    Nx, p = size(u)
    u_short = zeros(sL + sR, p)
    i = 1

    for k in j-sL+1:j+sR
        u_short[i,:] = u[mod1(k, Nx),:]
        i += 1
    end

    return u_short
end

function extractExtendedLocalData(u, j, sL, sR)

    # To keep only cells for which the modified flux is not G(K)
    Nx, p = size(u)

    u_short = zeros(2*(sR+sL),p)
    i=1
    for k in j-sL-sR+1:j+sR+sL
        u_short[i,:] = u[mod1(k, Nx), :]
        i+=1
    end

    return u_short
end

extractExtendedLocalData(u::Nothing, j, sL, sR) = nothing

computeK(clModifiedData::CLModifiedData, u) = sum(clModifiedData.weights .* u) / length(u)
function computeK(clModifiedData::CLModifiedData, u)
    Nx, p = size(u)
    K = zeros(p)
    for k in 1:p
        K[k] = sum(clModifiedData.weights .* u[:,k]) / Nx
    end
    K
end

#computeK(::MaxModifiedData, u) = max(u...)
computeK(::MaxModifiedData, u) = transpose(maximum(u, dims=1))
computeK(::MinModifiedData, u) = min(u...)

meanK(sL, sR) = CLModifiedData(ones(sL + sR))
maxK() = MaxModifiedData()
minK() = MinModifiedData()

computeZ(KFun::ModifiedDataType, z, j, sL, sR) = computeK(KFun, extractLocalData(z, j, sL, sR))
computeZ(::ModifiedDataType, ::Nothing, j, sL, sR) = nothing

abstract type BoundsType end
struct NormalBounds <: BoundsType end
struct LightBounds <: BoundsType end

function compute_u_tilde(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int)
    
    K = computeK(KFun, extractLocalData(u, j, sL, sR))
    #@show size(K)
    Nx, p = size(u) 

    ut = zero(u)

    #ut = K * ones(Nx, p)

    
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

    Nx = length(u)

    ut = zeros(Nx)

    ileft, iright = j - sL + 1, j + sR
    jleft, jright = j - sL - sR + 2, j + sL + sR - 1

    for k in 1:Nx

        if (
            (ileft <= k && iright >= k)
            || (ileft % Nx <= k && ileft < 1)
            || (iright % Nx >= k && iright > Nx)
        )
            ut[k] = u[k]
        elseif (jleft <= k && ileft > k) || (
            jleft % Nx <= k && jleft < 1 && k < ileft % Nx
        )
            ut[k] = u[mod1(ileft, Nx)]
        elseif (iright < k && jright >= k) || (
            jright % Nx >= k && jright > Nx && k > iright % Nx
        )
            ut[k] = u[mod1(iright, Nx)]
        elseif jright < k
            ut[k] = u[mod1(iright, Nx)]
        elseif jleft > k
            ut[k] = u[mod1(ileft, Nx)]
        end
    end
    ut
end

compute_z_tilde(::NullSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR) = nothing
function compute_z_tilde(zbSource::ZbSource, modifiedDataType::ModifiedDataType, domain::Domain, j, sL, sR)
    z = zeros(domain.Nx,1)
    for k in 1:Nx
        z[k] = zb(zbSource, domain.x[k])
    end
    compute_u_tilde(modifiedDataType, z, j, sL, sR)
end

function compute_u_hat(::NullSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=nothing)

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    for k in j-sL-sR+1:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            numFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :])
            .-
            numFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :])
        )

    end

    return uh
end

function compute_u_hat(::ZbSource, ut, dx, dt, j, domain::Domain, equation::Equation, method::FVMethod; zt=zero(ut))

    uh = copy(ut)
    Nx, p = size(ut)

    sL, sR = get_sL(method), get_sR(method)

    sourceVec = sourceTerm(method, domain, ut; z=zt)

    for k in j-sL-sR+1:j+sR+sL

        uh[mod1(k, Nx), :] = ut[mod1(k, Nx), :] .- dt / dx .* (
            giveNumFlux(method, equation, ut[mod1(k, Nx), :], ut[mod1(k + 1, Nx), :]; zL=zt[mod1(k, Nx)], zR=zt[mod1(k + 1, Nx)])
            .-
            giveNumFlux(method, equation, ut[mod1(k - 1, Nx), :], ut[mod1(k, Nx), :]; zL=zt[mod1(k - 1, Nx)], zR=zt[mod1(k, Nx)])) .+ dt * sourceVec[mod1(k, Nx),:]
        

    end

    return uh
end

# function compute_u_hat(ut, j::Int, domain::Domain, equation::Equation, method::FVMethod, zt=nothing)
#     uh = copy(ut)
#     Nx, p = size(ut)
#     sL, sR = get_sL(method), get_sR(method)

#     uh_short = extractExtendedLocalData(uh, j, sL, sR)
#     zt_short = extractExtendedLocalData(zt, j, sL, sR)


    
# end

function initBounds(KFun::SymmetricModifiedData, equation::Equation, u, j, sL, sR, z=nothing)
    GK = get_G(equation, computeK(KFun, extractLocalData(u, j, sL, sR)); z=computeZ(KFun, z, j, sL, sR))
    return GK[1], GK[1]
end

initBounds(::AsymmetricModifiedData, equation::Equation, u, j, sL, sR) = G(equation, u[mod1(j + sR, length(u))]), G(equation, u[mod1(j - sL + 1, length(u))])

function updateBounds!(::SymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt=zero(ut))

    for k in j-sL-sR+1:j
        M = M .+ dx / dt .* (get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]))
    end

    for k in j+1:j+sL+sR
        m = m .+ dx / dt .* (get_eta(equation, uh[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]) .- get_eta(equation, ut[mod1(k, Nx),:]; z=zt[mod1(k, Nx),:]))
        
    end
    
    m, M

end

function updateBounds!(::AsymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt)

    for k in j-sL-sR+2:j
        M += dx / dt * (eta(equation, ut[mod1(k, Nx)]) - eta(equation, uh[mod1(k, Nx)]))
    end

    for k in j+1:j+sL+sR-1
        m += dx / dt * (eta(equation, uh[mod1(k, Nx)]) - eta(equation, ut[mod1(k, Nx)]))
    end

    
    m, M

end

function compute_G_bounds(u, Nx, dx, dt, equation::Equation, domain::Domain, method::FVMethod, modifiedDataType::ModifiedDataType=SymmetricModifiedData(), boundsType::BoundsType=NormalBounds())

    M_vec, m_vec = zeros(Nx + 1), zeros(Nx + 1)
    sL, sR = get_sL(method), get_sR(method)

    # source
    # z = zeros(domain.Nx,1)
    # for i in eachindex(domain.x) z[i]=zb(equation.source, domain.x[i]) end
    z = isnothing(domain.sourceVec) ? zeros((Nx, 1)) : reshape(domain.sourceVec, (domain.Nx,1))

    for j in 1:Nx
        ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        # zt = compute_z_tilde(equation.source, modifiedDataType, domain, j, sL, sR)
        # zt = isnothing(zt) ? zero(ut) : zt
        zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
        uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)

        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR, z)

        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt, zt)

        # if m[1] > M[1]
        #     @warn "m greater than M !!!"
        # end

        m_vec[j+1], M_vec[j+1] = m[1], M[1]

    end

    m_vec[1], M_vec[1] = m_vec[end], M_vec[end]

    m_vec, M_vec
end

function J(gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation, domain::Domain)

    JD = 0
    JC = 0

    z = zb(equation.source, domain.x)

    for j in 1:Nx
        JD += max(0, get_eta(equation, up[j,:]; z=z[j])[1] - get_eta(equation, u[j,:]; z=z[j])[1] + dt / dx * (gamma[j+1] - gamma[j]))^2
    end

    for j in 1:Nx+1
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end

    JD + JC
end

function diffusion(u, up, gamma, dx, dt, equation, domain::Domain)
    z = zb(equation.source, domain.x)
    [get_eta(equation, up[i,:]; z=z[i])[1] - get_eta(equation, u[i,:]; z=z[i])[1] for i in 1:length(u[:,1])] + dt / dx * (gamma[2:end] - gamma[1:end-1])
end

initial_guess(m_vec, M_vec) = 0.5 * (m_vec + M_vec)

function optimize_for_entropy(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1, 1), boundsType::BoundsType=NormalBounds(), kwargs...)

    Nx, dx = domain.Nx, domain.dx
    FVsol = fv_solve(domain, u_init, equation, method)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, domain, method, modifiedDataType, boundsType)
    gamma_init = initial_guess(m_vec, M_vec)

    sol = optimize(gamma -> J(gamma, u_approx[end-1], u_approx[end], Nx, dx, dt_vec[end], m_vec, M_vec, equation, domain), gamma_init; iterations=100000, kwargs...)

    Gopt, Jopt = Optim.minimizer(sol), Optim.minimum(sol)
    Dopt = diffusion(u_approx[end-1], u_approx[end], Gopt, dx, dt_vec[end], equation, domain)
    OptForEntropySol(domain, equation, method, u_approx, dt_vec, Gopt, Jopt, Dopt, m_vec, M_vec)

end