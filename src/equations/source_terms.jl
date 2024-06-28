#include("abstract_equations.jl")
#include("../methods/rusanov.jl")
#include("../domain.jl")
#include("../methods/hydrostatic.jl")

abstract type Source end
struct NullSource <: Source end

abstract type ZbSource <: Source end
struct Bump_zb <: ZbSource end

zb(::Bump_zb, x) = (-0.5 .* x .* (-1 .+ x)) .* 5
Dzb(::Bump_zb, x) = (-x .+ 0.5) .* 5

sourceTerm(::FVMethod, zbSource::ZbSource, domain::Domain, v) = [[0.0, -v[i][1] * g * Dzb(zbSource, domain.x)[i]] for i in eachindex(v)]
