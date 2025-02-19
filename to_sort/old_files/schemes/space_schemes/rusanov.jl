struct Rusanov{T <: AbstractFloat} <: SpaceScheme
    CFL_factor::T
end

get_sL(::Rusanov) = 1
get_sR(::Rusanov) = 1

get_name(::Rusanov) = "Rusanov"

function numFlux(::Rusanov, equation::Equation, u; kwargs...)
    #@show [uL, uR]
    uL, uR = u[1,:], u[2,:]
    A = CFL_cond(equation, u)
    #@show (flux(equation, uL) .+ flux(equation, uR)) / 2 .- A / 2 * (uR .- uL)
    integrator.flux .= (flux(equation, uL) .+ flux(equation, uR)) / 2 .- A / 2 * (uR .- uL)
end

# function giveNumFlux(::Rusanov, equation::Equation, v; kwargs...)
#     Nx, p = length(v[:,1]), get_unknowns_number(equation)
#     numericalFluxMat = zeros(Nx+1, p)
#     for i ∈ 2:Nx
#         numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:]; zL=domain.sourceVec[i-1], zR=domain.sourceVec[i])
#     end
#     numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:]; zL=domain.sourceVec[end], zR=domain.sourceVec[1])
#     numericalFluxMat[end,:] = numericalFluxMat[1,:]
#     numericalFluxMat

# end

function exactEntropicNumFlux(::Rusanov, equation::Equation, uL, uR)
    A = max(abs(D_flux(equation, uL)), abs(D_flux(equation, uR)))
    return (G(equation, uL) + G(equation, uR)) / 2 - A / 2 * (eta(equation, uR) - eta(equation, uL))
end

function exactG(scheme::Rusanov, equation::Equation, u)
    Nx = length(u)
    G = zeros(Nx + 1)
    for i in 1:Nx-1
        G[i+1] = exactEntropicNumFlux(scheme, equation, u[i], u[i+1])
    end
    G[begin] = exactEntropicNumFlux(scheme, equation, u[end], u[begin])
    G[end] = G[begin]
    G
end