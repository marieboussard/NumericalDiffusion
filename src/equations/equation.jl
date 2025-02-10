# struct Equation <: AbstractEquation
#     p::Int
#     flux
#     Dflux
#     initcond
# end

struct Equation{T<:EquationType, equationFunType<:AbstractEquationFun, F1<:Base.Callable, sourceType<:AbstractSource} <: AbstractEquation{T}
    p::Int
    eqtype::T
    funcs::equationFunType
    initcond::F1
    source::sourceType
end

function Equation(p::Int, eqtype::T, funcs::equationFunType, initcond::F1) where T<: EquationType where equationFunType<:AbstractEquationFun where F1<:Base.Callable
    source = NoSource()
    Equation(p, eqtype, funcs, initcond, source)
    #{typeof(eqtype), typeof(funcs), typeof(source), typeof(initcond)}
end