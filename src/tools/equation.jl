#include("method.jl")

abstract type Equation end
abstract type ScalarEquation <: Equation end

# Handling of source terms

abstract type Source end
struct NullSource <: Source end



struct Burgers <: ScalarEquation
    source::Source
end

burgers() = Burgers(NullSource())

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

""" Saint Venant equation """

struct SaintVenant <: Equation
    source::Source
end

abstract type ZbSource <: Source end
struct Bump_zb <: ZbSource end

g = 9.8


function flux(::SaintVenant, v)
    h, hu = v
    return [hu, hu^2 / h + g * h^2 / 2]
end


function D_flux(::SaintVenant, v)

    res = Matrix[]
    for vi in v
        h, hu = vi
        push!(res, [0 1; -(hu^2)/h^2+g*h 2*hu/h])
        #res.append(np.array([[0, 1], [-(hu^2) / h^2 + g * h, 2 * hu / h]]))
    end
    return res
end



#bump_source(v, x) = [[0.0, -v[i][0] * g * (-x[i] + 0.5)] for i in eachindex(v)]

CFL_cond(equation::ScalarEquation, v) = max(abs.(D_flux(equation, v))...)
CFL_cond(::SaintVenant, v) = max([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]...)

# function CFL_cond(::SaintVenant, v)
#     @show typeof(v)
#     @show [vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]
#     @show typeof([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v])
#     @show max([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]...)
# end