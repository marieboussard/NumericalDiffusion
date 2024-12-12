#include("abstract_methods.jl")

struct Roe{T <: AbstractFloat} <: SpaceScheme
    CFL_factor::T
end
get_sL(::Roe) = 1
get_sR(::Roe) = 1

get_name(::Roe) = "Roe"

function numFlux(::Roe, equation::ScalarEquation, uL, uR)
    sigma = (flux(equation, uR) .- flux(equation, uL)) ./ (uR .- uL)
    if sigma[1] < 0
        return flux(equation, uR)
    else
        return flux(equation, uL)
    end
end