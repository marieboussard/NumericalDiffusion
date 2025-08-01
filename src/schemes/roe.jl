mutable struct RoeCache{utype<:AbstractArray} <: SpaceCache
    sigma::Float64
    fcont::utype
    function RoeCache(u::AbstractArray)
        fcont = zero(u)
        new{typeof(u)}(zero(Float64), fcont)
    end
end

struct Roe <: SpaceScheme end

get_sL(::Roe) = 1
get_sR(::Roe) = 1

get_name(::Roe) = "Roe"


function update_cache!(rcache::RoeCache, u::AbstractArray, equation::Equation, jstart::Int=1, jend::Int=length(u))
    @unpack fcont = rcache
    flux!(equation.funcs, view(u, jstart:jend), fcont)
end


function numflux!(::Roe, rcache::RoeCache, equation::Equation, u::AbstractVector, fnum::AbstractVector, ju::Int, Nx=length(u), jf::Int=ju)
    @unpack fcont = rcache

    for r in 1:equation.p
        rcache.sigma = (fcont[mod1(ju+1,Nx), r] - fcont[ju, r]) / (u[mod1(ju+1,Nx),r] - u[ju,r])
        if rcache.sigma<0
            fnum[jf,r] = fcont[mod1(ju+1,Nx), r]
        else
            fnum[jf,r] = fcont[ju, r]
        end
    end
end

# function numflux!(::Roe, integrator::Integrator, j::Int, args...)
#     @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
#     @unpack Nx = integrator.params.mesh
#     for r in 1:equation.p
#         space_cache.sigma = (fcont[mod1(j+1,Nx), r] - fcont[j, r]) / (uprev[mod1(j+1,Nx),r] - uprev[j,r])
#         if space_cache.sigma<0
#             fnum[j,r] = fcont[mod1(j+1,Nx), r]
#         else
#             fnum[j,r] = fcont[j, r]
#         end
#     end
# end

# function numflux!(::Roe, j::Int, params::Parameters, equation::Equation, cache::Cache, space_cache::SpaceCache, fnum::AbstractArray, fcont::AbstractArray, u::AbstractArray, i::Int=j)
#     @unpack Nx = params.mesh
#     for r in 1:equation.p
#         space_cache.sigma = (fcont[mod1(j+1,Nx), r] - fcont[j, r]) / (u[mod1(j+1,Nx),r] - u[j,r])
#         if space_cache.sigma<0
#             fnum[j,r] = fcont[mod1(j+1,Nx), r]
#         else
#             fnum[j,r] = fcont[j, r]
#         end
#     end
# end