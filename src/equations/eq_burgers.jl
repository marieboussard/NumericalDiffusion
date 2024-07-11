#include("abstract_equations.jl")

struct Burgers <: ScalarEquation
    source::Source
end

get_name(::Burgers) = "Burgers"

burgers() = Burgers(NullSource())

flux(::Burgers, u) = u.^2 ./ 2
D_flux(::Burgers, u) = u
eta(::Burgers, u) = u.^2
#etaTilde(::Burgers, u; x=0) = eta(::Burgers, u; x=0)
G(::Burgers, u) = 2.0 .* u.^3 ./ 3.0
#GTilde(::Burgers, u; x=0) = G(::Burgers, u; x=0)

get_eta(equation::Burgers, u; z=nothing) = eta(equation, u)
get_G(equation::Burgers, u; z=nothing) = G(equation, u)

# Initial conditions

function uexact_burgers_article(x::Real, t::Real)
    if t >= 2 / 3
        @warn "Warning: This solution is not valid for t ≥ 2/3"
    end
    if x <= -2 * t
        return -1 * (x + 2) / (1 - t)
    elseif x <= 3 * t
        return x / t
    else
        return -3 * (x - 2) / (2 - 3 * t)
    end
end

function u0_burgers_article(x::Real)
    if x <= 0
        return -2 - x
    else
        return 3 - 3 / 2 * x
    end
   
end
u0_burgers_article(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x