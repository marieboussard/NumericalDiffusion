abstract type Equation end
abstract type ScalarEquation <: Equation end
abstract type NoSourceScalarEquation <: ScalarEquation end

flux(equation::Equation, args...) = @error "Flux not given for equation "*get_name(equation)
eta(equation::Equation, args...) = @error "η not given for equation "*get_name(equation)
G(equation::Equation, args...) = @error "G not given for equation "*get_name(equation)
D_eta(equation::Equation, args...) = @error "Dη not given for equation "*get_name(equation)

# Functions for scalar equations
CFL_cond(equation::ScalarEquation, v::Vector{Vector{Real}}) = max(abs.(D_flux(equation, [v[i][1] for i in eachindex(v)]))...)
CFL_cond(equation::ScalarEquation, v::Vector{Vector{Float64}}) = max(abs.(D_flux(equation, [v[i][1] for i in eachindex(v)]))...)
#CFL_cond(equation::ScalarEquation, v::Vector{Float64}) = max(abs.(D_flux(equation, v))...)
CFL_cond(equation::ScalarEquation, v) = maximum(abs.(D_flux(equation, v)))
get_unknowns_number(::ScalarEquation) = 1

# ========== Testcases ========= #
abstract type Testcase end
u0_fun(testcase::Testcase, x) = "Missing function"

