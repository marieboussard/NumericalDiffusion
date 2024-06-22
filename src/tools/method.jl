include("equation.jl")

abstract type FVMethod end

struct Rusanov <: FVMethod
    CFL_factor::Real
end

struct Roe <: FVMethod
    CFL_factor::Real
end

function numFlux(::Rusanov, equation::Equation, uL::Real, uR::Real)

    A = max(abs(D_flux(equation, uL)), abs(D_flux(equation, uR)))
    (flux(equation, uL) + flux(equation, uR)) / 2 - A / 2 * (uR - uL)

end

function numFlux(::Roe, equation::Equation, uL::Real, uR::Real)

    A = max(abs(D_flux(equation, uL)), abs(D_flux(equation, uR)))
    (G(equation, uL) + G(equation, uR)) / 2 - A / 2 * (eta(equation, uR) - eta(equation, uL))

end