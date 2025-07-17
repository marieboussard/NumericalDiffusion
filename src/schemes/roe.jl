mutable struct RoeCache <: SpaceCache
    sigma::Float64
    RoeCache() = new(zero(Float64))
end

struct Roe <: SpaceScheme end

get_sL(::Roe) = 1
get_sR(::Roe) = 1

get_name(::Roe) = "Roe"

function numflux!(::Roe, integrator::Integrator, j::Int, args...)
    @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack Nx = integrator.params.mesh
    for r in 1:equation.p
        space_cache.sigma = (fcont[mod1(j+1,Nx), r] - fcont[j, r]) / (uprev[mod1(j+1,Nx),r] - uprev[j,r])
        if space_cache.sigma<0
            fnum[j,r] = fcont[mod1(j+1,Nx), r]
        else
            fnum[j,r] = fcont[j, r]
        end
    end
end

function numflux!(::Roe, j::Int, params::Parameters, equation::Equation, cache::Cache, space_cache::SpaceCache, fnum::AbstractArray, fcont::AbstractArray, u::AbstractArray, i::Int=j)
    @unpack Nx = params.mesh
    for r in 1:equation.p
        space_cache.sigma = (fcont[mod1(j+1,Nx), r] - fcont[j, r]) / (u[mod1(j+1,Nx),r] - u[j,r])
        if space_cache.sigma<0
            fnum[j,r] = fcont[mod1(j+1,Nx), r]
        else
            fnum[j,r] = fcont[j, r]
        end
    end
end