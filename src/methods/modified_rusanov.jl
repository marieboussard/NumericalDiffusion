struct ModifiedRusanov <: FVMethod
    CFL_factor::Float64
    A::Real
end

get_sL(::ModifiedRusanov) = 1
get_sR(::ModifiedRusanov) = 1

get_name(mr::ModifiedRusanov) = "Modified Rusanov" # (A="*string(mr.A)*")"

function numFlux(mr::ModifiedRusanov, equation::Equation, uL, uR)
    #[uL, uR]
    A = mr.A
    (flux(equation, uL) .+ flux(equation, uR)) / 2 .- A / 2 * (uR .- uL)
end