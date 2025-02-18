struct Advection <: AbstractEquationFun 
a::Float64
end

flux(eq::Advection, u) = eq.a*u
Dflux(eq::Advection, u) = eq.a

u0_gauss(x::Real; xm=0.0, sigma=0.1) = 1/sqrt(2*pi*sigma)*exp(-(x-xm)^2/sigma^2)
u0_gauss(x::AbstractArray) = u0_gauss.(x)

AdvectionExample = Equation(OneD(), 1, Scalar(), Advection(2.0), u0_gauss)

# u0_article(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
# u0_article(x::AbstractVector) = u0_article.(x)

# BurgersArticle = Equation(OneD(), 1, Scalar(), Burgers(), u0_article)