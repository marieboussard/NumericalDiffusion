#include("abstract_methods.jl")

struct Roe <: FVMethod
    CFL_factor::Float64
end
get_sL(::Roe) = 1
get_sR(::Roe) = 1

function numFlux(::Roe, equation::Equation, uL, uR)
    sigma = (flux(equation, uR) .- flux(equation, uL)) ./ (uR .- uL)
    if sigma < 0
        return flux(equation, uR)
    else
        return flux(equation, uL)
    end
end