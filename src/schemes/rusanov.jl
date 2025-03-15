mutable struct RusanovCache <: SpaceCache
    cfl_loc::Float64
    RusanovCache() = new(zero(Float64))
end

struct Rusanov <: SpaceScheme end

get_sL(::Rusanov) = 1
get_sR(::Rusanov) = 1

function numflux(::Rusanov, equation, u, args...)
    @unpack flux = equation.funcs
    @views uL = u[1,:]
    @views uR = u[2,:]
    (flux(uL) .+ flux(uR)) ./ 2 - CFL_cond(u, equation)./ 2 * (uR .- uL)
end

function numflux!(::Rusanov, integrator::Integrator, j::Int, args...)
    @unpack equation, params, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack Nx = params.mesh
    CFL_local!(equation.dim, equation.eqtype, equation.funcs, integrator, j)
    for r in 1:equation.p
        # fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
        fnum[j,r] = (fcont[j, r] + fcont[mod1(j+1,Nx),r]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[mod1(j+1,Nx),r] - uprev[j,r])
    end
end

function numflux!(::Rusanov, j::Int, params::Parameters, equation::Equation, cache::Cache, space_cache::SpaceCache, fnum::AbstractArray, fcont::AbstractArray, u::AbstractArray, i::Int=j)
    @unpack Nx = params.mesh
    CFL_local!(equation.dim, equation.eqtype, equation.funcs, j, params, cache, space_cache)
    for r in 1:equation.p
        fnum[i,r] = (fcont[j, r] + fcont[mod1(j+1,Nx),r]) *0.5 - space_cache.cfl_loc.*0.5 * (u[mod1(j+1,Nx),r] - u[j,r])
    end
end

# function numflux!(::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, subcache::scacheType, integrator::Integrator, j::Int, args...)
#     @unpack equation, params, fnum = integrator
#     @unpack Nx = params.mesh
#     subcache.cfl_loc = max(abs(uL[2]/uL[1]) + sqrt(g*uL[1]), abs(uR[2]/uR[1]) + sqrt(g*uR[1]))
#     for r in 1:equation.p
#         fnum[j,r] = (fL[r] + fR[r]) *0.5 - subcache.cfl_loc.*0.5 * (uR[r] - uL[r])
#     end
# end

function numflux!(::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, fnum::AbstractVector, subcache::SpaceCache, equation::Equation, args...)
    # @unpack equation, params, fnum = integrator
    # @unpack Nx = params.mesh
    subcache.cfl_loc = max(abs(uL[2]/uL[1]) + sqrt(g*uL[1]), abs(uR[2]/uR[1]) + sqrt(g*uR[1]))
    for r in 1:equation.p
        fnum[r] = (fL[r] + fR[r]) *0.5 - subcache.cfl_loc.*0.5 * (uR[r] - uL[r])
    end
end

numflux!(scheme::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, subcache::SpaceCache, integrator::Integrator, j::Int, args...) = numflux!(scheme, uL, uR, fL, fR, view(integrator.fnum,j,:), subcache, integrator.equation, args...)