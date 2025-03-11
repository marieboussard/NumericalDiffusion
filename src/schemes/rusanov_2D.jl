mutable struct Rusanov2DCache <: scacheType
    cflx_loc::Float64
    cfly_loc::Float64
    Rusanov2DCache() = new(zero(Float64), zero(Float64))
end

struct Rusanov2D <: SpaceScheme end

get_sL(::Rusanov2D) = 1
get_sR(::Rusanov2D) = 1

function numflux!(::Rusanov2D, integrator::Integrator, j::Int, k::Int, args...)
    # Fill numerical fluxes. Fj+1/2,k = fnum[j,k]; Hj,k+1/2=hnum[j,k]
    @unpack equation, cache, space_cache, fnum, fcont, uprev, params = integrator
    # @unpack fnum, hnum = fnum 
    # @unpack fcont, hcont = fcont
    @unpack Nx = params.mesh
    CFL_local!(equation.dim, equation.eqtype, integrator, j, k)
    for r in 1:equation.p
        # fnum[j,k,r] = (fcont[j,k,r]+ fcont[mod1(j+1,Nx),k,r]) *0.5 - space_cache.cflx_loc.*0.5 * (uprev[mod1(j+1,Nx),k,r] - uprev[j,k,r])
        # hnum[j,k,r] = (hcont[j,k,r]+ hcont[j,mod1(k+1,Nx),r]) *0.5 - space_cache.cfly_loc.*0.5 * (uprev[j,mod1(k+1,Nx),r] - uprev[j,k,r])
        fnum[j,k,r,1] = (fcont[j,k,r,1]+ fcont[mod1(j+1,Nx),k,r,1]) *0.5 - space_cache.cflx_loc.*0.5 * (uprev[mod1(j+1,Nx),k,r] - uprev[j,k,r])
        fnum[j,k,r,2] = (fcont[j,k,r,2]+ fcont[j,mod1(k+1,Nx),r,2]) *0.5 - space_cache.cfly_loc.*0.5 * (uprev[j,mod1(k+1,Nx),r] - uprev[j,k,r])
    end
end