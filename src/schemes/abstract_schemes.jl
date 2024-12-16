#include("../equations/abstract_equations.jl")
abstract type Scheme end

#abstract type FVMethod<:Method end
abstract type SpaceScheme<:Scheme end
abstract type TimeScheme<:Scheme end

numFlux(::Scheme, ::Equation, vL, vR) = @error "Missing function"
# function giveNumFlux(scheme::Scheme, equation::Equation, vL, vR, args...; kwargs...)
#     # println("Entering simple giveNumFlux")
#     # println("type of vL : "*string(typeof(vL)))
#     numFlux(scheme, equation::Equation, vL, vR, args...; kwargs...)
# end

computeTimeStep(scheme::Scheme, equation::Equation, domain::Domain, u) = method.CFL_factor * domain.dx / CFL_cond(equation, u)