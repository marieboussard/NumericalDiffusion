include("equation.jl")

abstract type FVMethod end

struct Rusanov <: FVMethod
    CFL_factor::Float64
end

get_sL(::Rusanov) = 1
get_sR(::Rusanov) = 1

struct Roe <: FVMethod
    CFL_factor::Float64
end
get_sL(::Roe) = 1
get_sR(::Roe) = 1

function numFlux(::Rusanov, equation::Equation, uL::Real, uR::Real)
    A = max(abs(D_flux(equation, uL)), abs(D_flux(equation, uR)))
    (flux(equation, uL) + flux(equation, uR)) / 2 - A / 2 * (uR - uL)
end

function numFlux(::Roe, equation::Equation, uL::Real, uR::Real)
    sigma = (flux(equation, uR) - flux(equation, uL)) / (uR - uL)
    if sigma < 0
        return flux(equation, uR)
    else
        return flux(equation, uL)
    end
end


function exactEntropicNumFlux(::Rusanov, equation::Equation, uL, uR)
    A = max(abs(D_flux(equation, uL)), abs(D_flux(equation, uR)))
    return (G(equation, uL) + G(equation, uR)) / 2 - A / 2 * (eta(equation, uR) - eta(equation, uL))
end

function exactG(method::Rusanov, equation::Equation, u)
    Nx = length(u)
    G = zeros(Nx + 1)
    for i in 1:Nx-1
        G[i+1] = exactEntropicNumFlux(method, equation, u[i], u[i+1])
    end
    G[1] = exactEntropicNumFlux(method, equation, u[end], u[1])
    G[end] = G[1]
    G
end