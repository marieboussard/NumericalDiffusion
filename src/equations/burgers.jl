struct Burgers <: AbstractEquationFun end

flux(::Burgers, u) = 0.5*u.^2
Dflux(::Burgers, u) = u

u0_article(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
u0_article(x::AbstractVector) = u0_article.(x)

BurgersArticle = Equation(1, Scalar(), Burgers(), u0_article)