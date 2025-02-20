# 1D ADVECTION

struct Advection <: AbstractEquationFun 
a::Float64
end

flux(eq::Advection, u) = eq.a*u
Dflux(eq::Advection, u) = eq.a

u0_gauss(x::Real; xm=0.0, sigma=0.1) = 1/sqrt(2*pi*sigma)*exp(-(x-xm)^2/sigma^2)
u0_gauss(x::AbstractArray) = u0_gauss.(x)

AdvectionExample = Equation(OneD(), 1, Scalar(), Advection(2.0), u0_gauss)

# 2D ADVECTION

struct Advection2D <: AbstractEquationFun
    a::Float64
    b::Float64
end

flux_f(eq::Advection2D, u) = eq.a*u
flux_h(eq::Advection2D, u) = eq.b*u
Dflux_f(eq::Advection2D, u) = eq.a
Dflux_h(eq::Advection2D, u) = eq.b

u0_gauss2(x::Real, y::Real; xm=0.0, ym=0.0, sigmax=0.2, sigmay=0.2, A=2.0) = A*exp(-(x-xm)^2/(2*sigmax^2)-(y-ym)^2/(2*sigmay^2))

Advection2Example = Equation(TwoD(), 1, Scalar(), Advection2D(2.0, 1.0), u0_gauss2)