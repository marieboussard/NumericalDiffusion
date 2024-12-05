#include("../equations/abstract_equations.jl")
abstract type Method end

abstract type FVMethod<:Method end

numFlux(::Method, ::Equation, vL, vR) = @error "Missing function"
function giveNumFlux(method::Method, equation::Equation, vL, vR; kwargs...)
    # println("Entering simple giveNumFlux")
    # println("type of vL : "*string(typeof(vL)))
    numFlux(method, equation::Equation, vL, vR)
end

computeTimeStep(method::FVMethod, equation::Equation, domain::Domain, u) = method.CFL_factor * domain.dx / CFL_cond(equation, u)