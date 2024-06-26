abstract type Equation end
abstract type ScalarEquation <: Equation end

CFL_cond(equation::ScalarEquation, v) = max(abs.(D_flux(equation, v))...)
get_unknowns_number(::ScalarEquation) = 1
