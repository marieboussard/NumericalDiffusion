include("abstract_equations.jl")
include("../methods/rusanov.jl")
include("../methods/hydrostatic.jl")

abstract type Source end
struct NullSource <: Source end

abstract type ZbSource <: Source end
struct Bump_zb <: ZbSource end

zb(::Bump_zb, x) = (-0.5 .* x .* (-1 .+ x)) .* 5
Dzb(::Bump_zb, x) = (-x .+ 0.5) .* 5

sourceTerm(::Rusanov, zbSource::ZbSource, domain::Domain, v) = [[0.0, -v[i][1] * g * Dzb(zbSource, domain.x)[i]] for i in eachindex(v)]
function sourceTerm(::Hydrostatic, zbSource::ZbSource, domain::Domain, v)
    x, dx, Nx = domain.x, domain.dx, domain.Nx
    1 / dx .* [[0.0, g / 2 * ((hMinus(v[nmod(i, Nx)], x[nmod(i, Nx)], x[nmod(i + 1, Nx)], zbSource))^2 - (hPlus(v[nmod(i, Nx)], x[nmod(i - 1, Nx)], x[nmod(i, Nx)], zbSource))^2)] for i in 1:Nx]
end