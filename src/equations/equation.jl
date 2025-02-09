# struct Equation <: AbstractEquation
#     p::Int
#     flux
#     Dflux
#     initcond
# end

struct Equation{T<:EquationType, F1} <: AbstractEquation{T}
    p::Int
    eqtype::T
    funcs::EquationFun
    initcond::F1
end