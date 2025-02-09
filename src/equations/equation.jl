# struct Equation <: AbstractEquation
#     p::Int
#     flux
#     Dflux
#     initcond
# end

struct Equation{T<:EquationType, equationFunType<:AbstractEquationFun, F1} <: AbstractEquation{T}
    p::Int
    eqtype::T
    funcs::equationFunType
    initcond::F1
end