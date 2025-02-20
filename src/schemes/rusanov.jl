mutable struct RusanovCache <: scacheType
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
    CFL_local!(equation.eqtype, integrator, j)
    for r in 1:equation.p
        # fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
        fnum[j,r] = (fcont[j, r] + fcont[mod1(j+1,Nx),r]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[mod1(j+1,Nx),r] - uprev[j,r])
    end
end