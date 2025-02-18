mutable struct RoeCache <: scacheType
    sigma::Float64
    RoeCache() = new(zero(Float64))
end

struct Roe <: SpaceScheme end

compute_sL(::Roe) = 1
compute_sR(::Roe) = 1

function numflux!(::Roe, integrator::Integrator, i, args...)
    @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack stencil = cache

    for j in 1:equation.p
        space_cache.sigma = (fcont[stencil[2], j] - fcont[stencil[1], j]) / (uprev[stencil[2],j] - uprev[stencil[1],j])
        if space_cache.sigma<0
            fnum[i,j] = fcont[stencil[2], j]
        else
            fnum[i,j] = fcont[stencil[1], j]
        end
    end
end