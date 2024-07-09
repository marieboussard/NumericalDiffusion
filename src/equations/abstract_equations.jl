abstract type Equation end
abstract type ScalarEquation <: Equation end

function CFL_cond(equation::ScalarEquation, v::Vector{Vector{Real}})
    #@show typeof(v)
    #@show D_flux(equation, [v[i][1] for i in eachindex(v)])
    max(abs.(D_flux(equation, [v[i][1] for i in eachindex(v)]))...)
end
function CFL_cond(equation::ScalarEquation, v::Vector{Vector{Float64}})
    #@show typeof(v)
    #@show D_flux(equation, [v[i][1] for i in eachindex(v)])
    max(abs.(D_flux(equation, [v[i][1] for i in eachindex(v)]))...)
end
#CFL_cond(equation::ScalarEquation, v::Vector{Float64}) = max(abs.(D_flux(equation, v))...)
CFL_cond(equation::ScalarEquation, v) = max(abs.(D_flux(equation, v))...)
get_unknowns_number(::ScalarEquation) = 1
