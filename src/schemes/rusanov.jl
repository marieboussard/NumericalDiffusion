mutable struct RusanovCache <: scacheType
    # absDfcont
    # function RusanovCache(Nx::Int)
    #     new(zeros(Float64, Nx))
    # end
    cfl_loc::Float64
    RusanovCache() = new(zero(Float64))
end

struct Rusanov <: SpaceScheme end

compute_sL(::Rusanov) = 1
compute_sR(::Rusanov) = 1

function numflux(::Rusanov, equation, u, args...)
    @unpack flux = equation.funcs
    @views uL = u[1,:]
    @views uR = u[2,:]
    (flux(uL) .+ flux(uR)) ./ 2 - CFL_cond(u, equation)./ 2 * (uR .- uL)
end

function numflux!(::Rusanov, integrator::Integrator, i, args...)
    @unpack equation, cache, space_cache, fnum, fcont, uprev = integrator
    @unpack stencil = cache
    CFL_local!(equation.eqtype, integrator)
    for j in 1:equation.p
        # cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
        fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
    end
end