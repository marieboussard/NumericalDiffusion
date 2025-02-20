mutable struct RoeCache <: scacheType
    sigma::Float64
    RoeCache() = new(zero(Float64))
end

struct Roe <: SpaceScheme end

get_sL(::Roe) = 1
get_sR(::Roe) = 1

function numflux!(::Roe, integrator::Integrator, j::Int, args...)
    @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack Nx = integrator.params.mesh
    # @unpack stencil = cache

    for r in 1:equation.p
        space_cache.sigma = (fcont[mod1(j+1,Nx), r] - fcont[j, r]) / (uprev[mod1(j+1,Nx),r] - uprev[j,r])
        if space_cache.sigma<0
            fnum[j,r] = fcont[mod1(j+1,Nx), r]
        else
            fnum[j,r] = fcont[j, r]
        end
    end

    # for j in 1:equation.p
    #     space_cache.sigma = (fcont[stencil[2], j] - fcont[stencil[1], j]) / (uprev[stencil[2],j] - uprev[stencil[1],j])
    #     if space_cache.sigma<0
    #         fnum[i,j] = fcont[stencil[2], j]
    #     else
    #         fnum[i,j] = fcont[stencil[1], j]
    #     end
    # end
end