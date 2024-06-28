#include("abstract_methods.jl")

struct Rusanov <: FVMethod
    CFL_factor::Float64
end

get_sL(::Rusanov) = 1
get_sR(::Rusanov) = 1

function numFlux(::Rusanov, equation::Equation, uL, uR)
    A = CFL_cond(equation, [uL, uR])
    (flux(equation, uL) .+ flux(equation, uR)) / 2 .- A / 2 * (uR .- uL)
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
    G[begin] = exactEntropicNumFlux(method, equation, u[end], u[begin])
    G[end] = G[begin]
    G
end