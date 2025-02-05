struct NewEq <: NoSourceScalarEquation
    source::Source
end

get_name(::NewEq) = "New"

newEq() = NewEq(NullSource())

flux(::NewEq, u) = u.^3 .+ u
D_flux(::NewEq, u) = 3*u.^2 .+ 1
eta(::NewEq, u) = u.^2
G(::NewEq, u) = u.^2 .+ (3.0./2.0.*u.^4)
D_eta(::NewEq, u) = 2.0*u


# Initial conditions

# function uexact_burgers_article(x::Real, t::Real)
#     if t >= 2 / 3
#         @warn "Warning: This solution is not valid for t ≥ 2/3"
#     end
#     if x <= -2 * t
#         return -1 * (x + 2) / (1 - t)
#     elseif x <= 3 * t
#         return x / t
#     else
#         return -3 * (x - 2) / (2 - 3 * t)
#     end
# end

# function u0_con(x::Real)
#     if x <= 0
#         return 1
#     else
#         return -0.75
#     end
   
# end
u0_concave_convex(x::Real) = x <= 0 ? 1.0 : -0.75

function uexact_concave_convex(x::Real, t::Real)
    ug, ud = 1, -0.75
    s = 1 + ud^2 + ug^2 + ud*ug
    x ≤ t*s ? ug : ud
end

# Testcase of a non classical shock by a concave-convex flux

struct ConcaveConvexTestcase <: Testcase
end

u0_fun(::ConcaveConvexTestcase, x::Real) = u0_concave_convex(x)
uexact_fun(::ConcaveConvexTestcase, x::Real, t::Real) = uexact_concave_convex(x, t)