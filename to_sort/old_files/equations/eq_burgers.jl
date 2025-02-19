struct Burgers <: NoSourceScalarEquation
    source::Source
end

get_name(::Burgers) = "Burgers"

burgers() = Burgers(NullSource())

flux(::Burgers, u, args...; kwargs...) = u.^2 ./ 2
D_flux(::Burgers, u, args...; kwargs...) = u
eta(::Burgers, u, args...; kwargs...) = u.^2
G(::Burgers, u, args...; kwargs...) = 2.0 .* u.^3 ./ 3.0
D_eta(::Burgers, u, args...; kwargs...) = 2.0*u

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

# function u0_burgers_article(x::Real)
#     if x <= 0
#         return -2 - x
#     else
#         return 3 - 3 / 2 * x
#     end
   
# end
u0_burgers_article(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x

u0_simple(x::Real) = x<=0 ? -1 : 1

# Testcase for article: A posteriori entropy inequality

struct ArticleTestcase <: Testcase
end

u0_fun(::ArticleTestcase, x::Real) = u0_burgers_article(x)
uexact_fun(::ArticleTestcase, x::Real, t::Real) = uexact_burgers_article(x, t)

# A simple testcase to study properties of multidimensionnal a priori estimation of the diffusion
struct SimpleShock <: Testcase end
u0_fun(::SimpleShock, x::Real) = u0_simple(x)

# A testcase for a counter example
struct CounterExample <: Testcase end
u0_counter(x::Real) = x<=0 ? -2-3*x : 7 - 3/2*x
function u_exact_counter(x::Real, t::Real)
    #if t >= 2 / 7
    if t >= 14/9
        @warn "Warning: This solution is not valid for t ≥ 14/9"
    end
    if x <= -2*t
        return (3*x+2) / (3*t-1)
    elseif x >= 7*t
        return (3*x-14) / (3*t-2)
    else
        return x/t
    end
end
u0_fun(::CounterExample, x::Real) = u0_counter(x)
uexact_fun(::CounterExample, x::Real, t::Real) = u_exact_counter(x,t)

# A smooth testcase
struct SmoothTestcase<: Testcase end
u0_fun(::SmoothTestcase, x::Real) = uexact_burgers_article(x, 0.1)
uexact_fun(::SmoothTestcase, x::Real, t::Real) = uexact_burgers_article(x, t+0.1)

# A parametrized piecewise linear testcase.
struct PiecewiseLinear<: Testcase
    a1
    b1 
    a2 
    b2
end
u0_fun(pl::PiecewiseLinear, x::Real) =(@unpack a1, b1, a2, b2 = pl; x<=0 ? -a1*x-b1 : -a2*x+b2)
function uexact_fun(pl::PiecewiseLinear, x::Real, t::Real)
    @unpack a1, b1, a2, b2 = pl
    if t >= min(1/a1, 1/a2)
        @warn "Warning: This solution is not valid for t ≥ min(1/a1, 1/a2)"
    end
    if x <= -b1*t
        return (a1*x+b1) / (a1*t-1)
    elseif x >= b2*t
        return (a2*x-b2) / (a2*t-1)
    else
        return x/t
    end
end

spaceBounds(pl::PiecewiseLinear) = (@unpack a1, b1, a2, b2 = pl; (-b1/a1, b2/a2))
function integerNx(Nx::Int, pl::PiecewiseLinear)
    @unpack a1, b1, a2, b2 = pl
    Int(ceil((1 + floor(Nx/(a2*b1)))*(a2*b1 + a1*b2)*a2*b1))
end