using Optim

include("./tools/fv_solution.jl")
include("./tools/solution.jl")
#include("./tools/utils.jl")

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
    Nx = length(u)
    u_short = zeros(sL + sR)
    i = 1

    for k in j-sL+1:j+sR
        u_short[i] = u[nmod(k, Nx)]
        i += 1
    end

    return u_short
end

computeK(clModifiedData::CLModifiedData, u) = sum(clModifiedData.weights .* u) / length(u)
computeK(::MaxModifiedData, u) = max(u...)
computeK(::MinModifiedData, u) = min(u...)

meanK(sL, sR) = CLModifiedData(ones(sL + sR))
maxK() = MaxModifiedData()
minK() = MinModifiedData()

abstract type BoundsType end
struct NormalBounds <: BoundsType end
struct LightBounds <: BoundsType end

function compute_u_tilde(KFun::SymmetricModifiedData, u, j::Int, sL::Int, sR::Int; K=0)

    K = computeK(KFun, extractLocalData(u, j, sL, sR))

    Nx = length(u)

    ut = K * ones(Nx)

    for i in 1:Nx
        if (j - sL + 1 <= i && i <= j + sR) || (((j - sL + 1) % Nx <= i) && (j - sL + 1 < 1)) || ((i <= (j + sR) % Nx) & (j + sR > Nx))
            ut[i] = u[i]
        end
    end
    ut

end


function compute_u_tilde(::AsymmetricModifiedData, u, j::Int, sL::Int, sR::Int; K=0)

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
            ut[k] = u[nmod(ileft, Nx)]
        elseif (iright < k && jright >= k) || (
            jright % Nx >= k && jright > Nx && k > iright % Nx
        )
            ut[k] = u[nmod(iright, Nx)]
        elseif jright < k
            ut[k] = u[nmod(iright, Nx)]
        elseif jleft > k
            ut[k] = u[nmod(ileft, Nx)]
        end
    end
    ut
end

function compute_u_hat(ut, dx, dt, j, equation, method::FVMethod)

    uh = copy(ut)
    Nx = length(ut)

    sL, sR = get_sL(method), get_sR(method)

    for k in j-sL-sR+1:j+sR+sL

        uh[nmod(k, Nx)] = ut[nmod(k, Nx)] - dt / dx * (
            numFlux(method, equation, ut[nmod(k, Nx)], ut[nmod(k + 1, Nx)])
            -
            numFlux(method, equation, ut[nmod(k - 1, Nx)], ut[nmod(k, Nx)])
        )

    end

    return uh
end

function initBounds(KFun::SymmetricModifiedData, equation::Equation, u, j, sL, sR)
    GK = G(equation, computeK(KFun, extractLocalData(u, j, sL, sR)))
    return GK, GK
end

initBounds(::AsymmetricModifiedData, equation::Equation, u, j, sL, sR) = G(equation, u[nmod(j + sR, length(u))]), G(equation, u[nmod(j - sL + 1, length(u))])

function updateBounds!(::SymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt)

    for k in j-sL-sR+1:j
        M += dx / dt * (eta(equation, ut[nmod(k, Nx)]) - eta(equation, uh[nmod(k, Nx)]))
    end

    for k in j+1:j+sL+sR
        m += dx / dt * (eta(equation, uh[nmod(k, Nx)]) - eta(equation, ut[nmod(k, Nx)]))
    end

    m, M

end

function updateBounds!(::AsymmetricModifiedData, ::NormalBounds, equation::Equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt)

    for k in j-sL-sR+2:j
        M += dx / dt * (eta(equation, ut[nmod(k, Nx)]) - eta(equation, uh[nmod(k, Nx)]))
    end

    for k in j+1:j+sL+sR-1
        m += dx / dt * (eta(equation, uh[nmod(k, Nx)]) - eta(equation, ut[nmod(k, Nx)]))
    end

    m, M

end

function compute_G_bounds(u, Nx, dx, dt, equation::Equation, method::FVMethod, modifiedDataType::ModifiedDataType=SymmetricModifiedData(), boundsType::BoundsType=NormalBounds())

    M_vec, m_vec = zeros(Nx + 1), zeros(Nx + 1)
    sL, sR = get_sL(method), get_sR(method)

    for j in 1:Nx
        ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        uh = compute_u_hat(ut, dx, dt, j, equation, method)

        m, M = initBounds(modifiedDataType, equation, u, j, sL, sR)
        m, M = updateBounds!(modifiedDataType, boundsType, equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt)

        if m > M
            @warn "m greater than M !!!"
        end

        m_vec[j+1], M_vec[j+1] = m, M

    end

    m_vec[1], M_vec[1] = m_vec[end], M_vec[end]

    m_vec, M_vec
end

function J(gamma, u, up, Nx, dx, dt, m_vec, M_vec, equation::Equation)

    JD = 0
    JC = 0

    for j in 1:Nx
        JD += max(0, eta(equation, up[j]) - eta(equation, u[j]) + dt / dx * (gamma[j+1] - gamma[j]))^2
    end

    for j in 1:Nx+1
        JC += (dt / dx)^2 * (max(0, (gamma[j] - M_vec[j]))^2 + max(0, 1 * (m_vec[j] - gamma[j]))^2)
    end

    JD + JC
end

diffusion(u, up, gamma, dx, dt, equation) = [eta(equation, up[i]) - eta(equation, u[i]) for i in eachindex(u)] + dt / dx * (gamma[2:end] - gamma[1:end-1])

initial_guess(m_vec, M_vec) = 0.5 * (m_vec + M_vec)

function optimize_for_entropy(u_init, domain::Domain, equation::Equation, method::FVMethod; modifiedDataType::ModifiedDataType=meanK(1, 1), boundsType::BoundsType=NormalBounds(), kwargs...)

    Nx, dx = domain.Nx, domain.dx
    FVsol = fv_solve(domain, u_init, equation, method)
    u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec

    m_vec, M_vec = compute_G_bounds(u_approx[end-1], Nx, dx, dt_vec[end], equation, method, modifiedDataType, boundsType)
    gamma_init = initial_guess(m_vec, M_vec)

    sol = optimize(gamma -> J(gamma, u_approx[end-1], u_approx[end], Nx, dx, dt_vec[end], m_vec, M_vec, equation), gamma_init; iterations=100000, kwargs...)

    Gopt, Jopt = Optim.minimizer(sol), Optim.minimum(sol)
    Dopt = diffusion(u_approx[end-1], u_approx[end], Gopt, dx, dt_vec[end], equation)
    OptForEntropySol(domain, equation, method, u_approx, dt_vec, Gopt, Jopt, Dopt, m_vec, M_vec)

end