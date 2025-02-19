# For the moment, we deal only with the one dimensional advection equation
struct Advection <: NoSourceScalarEquation
    velocity::Real
    source::Source
end

get_name(::Advection) = "Advection"
advection(a=1.0) = Advection(a, NullSource())

flux(eq::Advection, u) = eq.velocity * u
D_flux(eq::Advection, u) = eq.velocity
eta(::Advection, u) = 0.5 * u.^2
G(eq::Advection, u) = 0.5 * eq.velocity .* u.^2
D_eta(::Advection, u) = u

# Initial condition

u0_gaussian(x) = exp.(-(x.-0.5).^2 .*400)