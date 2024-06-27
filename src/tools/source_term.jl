include("equation.jl")
include("method.jl")


zb(::Bump_zb, x) = (-0.5 .* x .* (-1 .+ x)) .* 5
Dzb(::Bump_zb, x) = (-x .+ 0.5) .* 5

# zb(::Bump_zb, x) = 0.1 * cos.(100 * pi * x)
# Dzb(::Bump_zb, x) = -0.1 * 100 * pi * sin.(10 * pi * x)

#perturb(x, xleft, xright) = sin(pi * (x - xleft) / (xright - xleft))
perturb(x, xleft, xright) = 1

function h_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.1)
    #@show xleft, xright, x
    if x < xleft
        return [max(0, c - zb(zbSource, x)), 0]
    elseif x <= xright
        return [max(0, c - zb(zbSource, x)) + c * height * perturb(x, xleft, xright), 0]
    else
        return [max(0, c - zb(zbSource, x)), 0]
    end
end


#v0_lake_at_rest(x, zbSource::ZbSource; c=1) = [c - zb(zbSource, x), 0]
v0_lake_at_rest(x, zbSource::ZbSource; c=1) = [[max(0, c - zb(zbSource, xi)), 0] for xi in x]

# function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, perturbationWide=0.1, perturbationHeight=0.2, perturbationLoc=0.25)
#     res = [[c - zb(zbSource, xi), 0] for xi in x]
#     height, semiWide, center = perturbationHeight * c, Int(round(perturbationWide * length(x) * 0.5)), Int(round(perturbationLoc * length(x)))

#     for k in center-semiWide:center
#         res[k][1] += height / semiWide * (k - (center - semiWide))
#     end

#     for k in center+1:center+semiWide
#         res[k][1] += -height / semiWide * (k - (center + semiWide))
#     end

#     res

# end

function v0_lake_at_rest_perturbated(x, zbSource::ZbSource; c=1, xleft=0.2, xright=0.3, height=0.4)
    return [h_perturbated(xi, zbSource; c=c, xleft=xleft, xright=xright, height=height) for xi in x]
end

sourceTerm(::Rusanov, zbSource::ZbSource, domain::Domain, v) = [[0.0, -v[i][1] * g * Dzb(zbSource, domain.x)[i]] for i in eachindex(v)]
function sourceTerm(::Hydrostatic, zbSource::ZbSource, domain::Domain, v)
    x, dx, Nx = domain.x, domain.dx, domain.Nx
    # hPlusMinus = hMinus(v[nmod(i, Nx)], x[nmod(i, Nx)], x[nmod(i + 1, Nx)], zbSource)
    # hMinusPlus = hPlus(v[nmod(i, Nx)], x[nmod(i - 1, Nx)], x[nmod(i, Nx)], zbSource)
    #hminus, hplus = hMinusPlus(vL, vR, xL, xR, zbSource::ZbSource)
    #dx .* [[0.0, g / 2 * ((hMinus(v[nmod(i, Nx)], x[nmod(i, Nx)], x[nmod(i + 1, Nx)], zbSource))^2 - (hPlus(v[nmod(i + 1, Nx)], x[nmod(i, Nx)], x[nmod(i + 1, Nx)], zbSource))^2)] for i in 1:Nx]
    1 / dx .* [[0.0, g / 2 * ((hMinus(v[nmod(i, Nx)], x[nmod(i, Nx)], x[nmod(i + 1, Nx)], zbSource))^2 - (hPlus(v[nmod(i, Nx)], x[nmod(i - 1, Nx)], x[nmod(i, Nx)], zbSource))^2)] for i in 1:Nx]
end