#include("../equations/abstract_equations.jl")

abstract type FVMethod end

numFlux(method::FVMethod, equation::Equation, vL, vR) = @error "Missing function"
giveNumFlux(method::FVMethod, equation::Equation, vL, vR; kwargs...) = numFlux(method, equation::Equation, vL, vR)
