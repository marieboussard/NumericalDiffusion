# include("abstract_equations.jl")
# include("source_terms.jl")

g = 9.8

struct SaintVenant <: Equation
    source::Source
end

flux(::SaintVenant, v) = (h = v[1]; hu = v[2]; [hu, hu^2 / h + g * h^2 / 2])

function D_flux(::SaintVenant, v)
    res = Matrix[]
    for vi in v
        h, hu = vi
        push!(res, [0 1; -(hu^2)/h^2+g*h 2*hu/h])
    end
    return res
end

eta(equation::SaintVenant, v, x) = v[2]^2/(2*v[1]) + g/2*v[1]^2 + v[1]*g*zb(equation.source, x)
G(equation::SaintVenant, v, x) = (v[2]^2/(2*v[1]) + g*v[1]^2 + v[1]*g*zb(equation.source, x))*v[2]/v[1]

CFL_cond(::SaintVenant, v::Vector) = max([vi[2] / vi[1] + sqrt(g * vi[1]) for vi in v]...)
CFL_cond(::SaintVenant, v::Matrix) = max((v[:,2]./v[:,1] + sqrt.(g*v[:,1]))...)

get_unknowns_number(::SaintVenant) = 2

"""Initial conditions"""

# 1 # Lake at rest 

#v0_lake_at_rest(x, zbSource::ZbSource; c=1) = [[max(0, c - zb(zbSource, xi)), 0] for xi in x]
function v0_lake_at_rest(x, zbSource::ZbSource; c=1)
    v0 = zeros(size(x)..., 2)
    for I in CartesianIndices(x)
        v0[I,1] = max(0, c - zb(zbSource, x[I]))
        v0[I,2] = 0
    end
    v0
end

# 2 # Lake at rest with a perturbation

perturb(x, xleft, xright) = 1
#perturb(x, xleft, xright) = sin(pi * (x - xleft) / (xright - xleft))

# function h_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.1)
#     if x < xleft
#         return [max(0, c - zb(zbSource, x)), 0]
#     elseif x <= xright
#         return [max(0, c - zb(zbSource, x)) + c * height * perturb(x, xleft, xright), 0]
#     else
#         return [max(0, c - zb(zbSource, x)), 0]
#     end
# end

function h_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.1)
    if x < xleft
        return max(0, c - zb(zbSource, x))
    elseif x <= xright
        return max(0, c - zb(zbSource, x)) + c * height * perturb(x, xleft, xright)
    else
        return max(0, c - zb(zbSource, x))
    end
end


# function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.4)
#     return [h_perturbated(xi, zbSource; c=c, xleft=xleft, xright=xright, height=height) for xi in x]
# end

function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.4)
    v0 = zeros(size(x)..., 2)
    for I in CartesianIndices(x)
        v0[I,1] = h_perturbated(x[I], zbSource; c=c, xleft=xleft, xright=xright, height=height)
        v0[I,2] = 0
    end
    v0
end

## TO DO ## maybe put Source instead of ZbSource