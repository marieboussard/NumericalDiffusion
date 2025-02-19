mutable struct Rusanov2DCache <: scacheType
end

struct Rusanov2D <: SpaceScheme end

function numflux!(::Rusanov2D, integrator::Integrator, i, args...)
    @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack stencil = cache
    CFL_local!(equation.eqtype, integrator)
    for j in 1:equation.p
        # cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
        fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
    end
end