#include("../equations/abstract_equations.jl")

abstract type FVMethod end

numFlux(method::FVMethod, equation::Equation, vL, vR) = @error "Missing function"
function giveNumFlux(method::FVMethod, equation::Equation, vL, vR; kwargs...)
    # println("Entering simple giveNumFlux")
    # println("type of vL : "*string(typeof(vL)))
    numFlux(method, equation::Equation, vL, vR)
end
