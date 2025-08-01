abstract type EquationDim end
struct OneD <: EquationDim end
struct TwoD <: EquationDim end

struct Equation{dimType<:EquationDim, T<:EquationType, equationFunType<:AbstractEquationFun, F1<:Base.Callable, sourceType<:AbstractSource} <: AbstractEquation{T}
    dim::dimType
    p::Int
    eqtype::T
    funcs::equationFunType
    initcond::F1
    source::sourceType
end

function Equation(dim::dimType, p::Int, eqtype::T, funcs::equationFunType, initcond::F1) where dimType<:EquationDim where T<: EquationType where equationFunType<:AbstractEquationFun where F1<:Base.Callable
    source = NoSource()
    Equation(dim, p, eqtype, funcs, initcond, source)
    #{typeof(eqtype), typeof(funcs), typeof(source), typeof(initcond)}
end