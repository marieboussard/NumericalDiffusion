abstract type Equation end
abstract type ScalarEquation <: Equation end

struct Burgers <: ScalarEquation end

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

struct SaintVenant <: Equation end

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


# def v0_lake_at_rest(x, zb, c):
#     # c is the constant surface (c=h+z)
#     return np.array([c - zb(x), np.zeros(len(x))]).T


# def v0_init(x):
#     return v0_lake_at_rest(x, lambda x: -0.5 * x * (x - 1), 1)


# def source_ex(v, x):
#     return np.array([np.zeros(v.T[0].shape), -v.T[0] * g * (-x + 0.5)]).T

abstract type ZbFun end
struct Bump_zb <: ZbFun end

zb(::Bump_zb, x) = -0.5 .* x .* (-1 .+ x)

#v0_lake_at_rest(x, zbFun::ZbFun; c=1) = [[c - zb(zbFun, xi), 0] for xi in x]
v0_lake_at_rest(x, zbFun::ZbFun; c=1) = [c - zb(zbFun, x), 0]

bump_source(v, x) = [[0.0, -v[i][0] * g * (-x[i] + 0.5)] for i in eachindex(v)]

CFL_cond(equation::ScalarEquation, v) = max(abs.(D_flux(equation, v))...)
CFL_cond(::SaintVenant, v) = max([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]...)

# function CFL_cond(::SaintVenant, v)
#     @show typeof(v)
#     @show [vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]
#     @show typeof([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v])
#     @show max([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]...)
# end