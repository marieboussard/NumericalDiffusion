abstract type AbstractEquationFun end

struct EquationFun{F1, F2} <: AbstractEquationFun
    flux::F1
    Dflux::F2
end

flux(eqfun::EquationFun, u) = eqfun.flux.(u)
Dflux(eqfun::EquationFun, u) = eqfun.Dflux.(u)