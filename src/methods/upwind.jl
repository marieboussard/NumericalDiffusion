struct Upwind <: FVMethod
    CFL_factor::Real
end

get_sL(::Upwind) = 1
get_sR(::Upwind) = 1

get_name(::Upwind) = "Upwind"

function numFlux(::Upwind, equation::Equation, uL, uR)
    A = 0.5*(Df(equation, uL) + Df(equation, uR))
    if A > 0
        return uL
    else
        return uR
    end
end