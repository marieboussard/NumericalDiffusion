# include("abstract_equations.jl")
# include("source_terms.jl")

g = 9.8

struct SaintVenant <: Equation
    source::Source
    treshold::Float64
end

function flux(equation::SaintVenant, v)
    h = v[1]
    if h > equation.treshold
        hu = v[2]; return [hu, hu^2 / h + g * h^2 / 2]
    else
        return [0.0, 0.0]
    end
end

function D_flux(equation::SaintVenant, v)
    res = Matrix[]
    for vi in v
        h, hu = vi
        if h > equation.treshold
            push!(res, [0 1; -(hu^2)/h^2+g*h 2*hu/h])
        else
            push!(res, [0 1; 0.0 0.0])
        end
    end
    return res
end

#eta(equation::SaintVenant, v, x) = v[2]^2/(2*v[1]) + g/2*v[1]^2 + v[1]*g*zb(equation.source, x)
function eta(equation::SaintVenant, v, z)
    if v[1] < equation.treshold
        return 0.0
    else
        return v[2]^2/(2*v[1]) .+ g/2*v[1]^2 .+ v[1]*g*z
    end
end
#etaTilde(equation::SaintVenant, v, x) = v[2]^2/(2*v[1]) + g/2*v[1]^2 + v[1]*g*zb(equation.source, x)
#G(equation::SaintVenant, v, x) = (v[2]^2/(2*v[1]) + g*v[1]^2 + v[1]*g*zb(equation.source, x))*v[2]/v[1]
function G(equation::SaintVenant, v, z)
    #@show v, z
    #@show (v[2]^2/(2*v[1]) .+ g*v[1]^2 .+ v[1]*g*z)*v[2]/v[1]

    if v[1] < equation.treshold
        return 0
    else
        return (v[2]^2/(2*v[1]) .+ g*v[1]^2 .+ v[1]*g.*z)*v[2]/v[1]
    end
end
get_eta(equation::SaintVenant, v; z=zero(v)) = eta(equation, v, z)
get_G(equation::SaintVenant, v; z=zero(v)) = G(equation, v, z)

function CFL_cond(equation::SaintVenant, v::Vector)
    n = size(v)[1]
    lamb = zeros(n)
    for i in 1:n
        if v[i][1] > equation.treshold
            lamb[i] = v[i][2] / v[i][1] + sqrt(g * abs.(v[i][1]))
        else
            lamb[i] = 0.0
        end
    end
    return maximum(lamb)
            
#CFL_cond(equation::SaintVenant, v::Vector) = max([vi[2] / max(vi[1], 1e-3) + sqrt(g * abs.(vi[1])) for vi in v]...)
end

function CFL_cond(equation::SaintVenant, v::Matrix)
    n = size(v)[1]
    lamb = zeros(n)
    for i in 1:n
        if v[i,1] > equation.treshold
            lamb[i] = v[i,2] / v[i,1] + sqrt(g * abs.(v[i,1]))
        else
            lamb[i] = 0.0
        end
    end
    return maximum(lamb)
end

#CFL_cond(::SaintVenant, v::Matrix) = max((v[:,2]./max.(v[:,1], 1e-3) + sqrt.(g*abs.(v[:,1])))...)

get_unknowns_number(::SaintVenant) = 2

"""Initial conditions"""

# 1 # Lake at rest 

#v0_lake_at_rest(x, zbSource::ZbSource; c=1) = [[max(0, c - zb(zbSource, xi)), 0] for xi in x]
function v0_lake_at_rest(x, source::Source; c=1)
    v0 = zeros(size(x)..., 2)
    for I in CartesianIndices(x)
        v0[I,1] = max(0, c - zb(source, x[I]))
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

    #@show c - zb(zbSource,x)
    if x < xleft
        #@show max(0, c - zb(zbSource, x))
        return max(0, c - zb(zbSource, x))
    elseif x <= xright
        #@show max(0, c - zb(zbSource, x))
        #@show c * height * perturb(x, xleft, xright)
        #@show max(0, c - zb(zbSource, x)) + c * height * perturb(x, xleft, xright) + zb(zbSource, x)
        return max(0, c - zb(zbSource, x)) + c * height * perturb(x, xleft, xright)
    else
        #@show max(0, c - zb(zbSource, x))
        return max(0, c - zb(zbSource, x))
    end
end


# function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.4)
#     return [h_perturbated(xi, zbSource; c=c, xleft=xleft, xright=xright, height=height) for xi in x]
# end

function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, xleft=0.1, xright=0.2, height=0.4)
    v0 = zeros(size(x)..., 2)
    for I in CartesianIndices(x)
        v0[I,1] = h_perturbated(x[I], zbSource; c=c, xleft=xleft, xright=xright, height=height)
        v0[I,2] = 0
    end
    v0
end

## TO DO ## maybe put Source instead of ZbSource