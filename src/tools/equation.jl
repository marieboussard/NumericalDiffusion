abstract type Equation end

struct Burgers <: Equation end

flux(::Burgers, x) = x^2 / 2
D_flux(::Burgers, x) = x
eta(::Burgers, x) = x^2
G(::Burgers, x) = 2.0 * x^3 / 3.0

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