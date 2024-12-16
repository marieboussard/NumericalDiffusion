struct Burgers <: NoSourceScalarEquation
    source::Source
end

get_name(::Burgers) = "Burgers"

burgers() = Burgers(NullSource())

flux(::Burgers, u) = u.^2 ./ 2
D_flux(::Burgers, u) = u
eta(::Burgers, u) = u.^2
G(::Burgers, u) = 2.0 .* u.^3 ./ 3.0
D_eta(::Burgers, u) = 2.0*u

# get_flux(equation::Burgers, u; z=nothing) = flux(equation, u)
# get_eta(equation::Burgers, u; z=nothing) = eta(equation, u)
# get_G(equation::Burgers, u; z=nothing) = G(equation, u)
# get_D_eta(equation::Burgers, u; z=nothing) = D_eta(equation, u)

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
        @warn "Warning: This solution is not valid for t ≥ 2/3"
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