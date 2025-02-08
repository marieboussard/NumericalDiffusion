struct RusanovCache <: scacheType
    # absDfcont
    # function RusanovCache(Nx::Int)
    #     new(zeros(Float64, Nx))
    # end
end

struct Rusanov <: SpaceScheme end

compute_sL(::Rusanov) = 1
compute_sR(::Rusanov) = 1

function numflux(::Rusanov, equation, u, args...)
    @unpack flux = equation
    @views uL = u[1,:]
    @views uR = u[2,:]
    (flux(uL) .+ flux(uR)) ./ 2 - CFL_cond(u, equation)./ 2 * (uR .- uL)
end

function numflux!(::Rusanov, integrator::Integrator, i, args...)
    @unpack equation, cache, fnum, fcont, uprev = integrator
    @unpack stencil = integrator.cache

    #=
    I = stencil[1]
    J = stencil[2]

    uL   = view(uprev, stencil[1])
    uR   = view(uprev, stencil[2])

    
    fL   = view(fcont, stencil[1])
    fR   = view(fcont, stencil[2])
    #fnum_i = view(fnum, i, :)
    =#
    cache.cfl_loc = maximum([abs.(uprev[stencil[1]]), abs.(uprev[stencil[2]])])
    fnum[i] = 0.5 * (fcont[stencil[1]] + fcont[stencil[2]]) - 0.5 * cache.cfl_loc * (uprev[stencil[2]] - uprev[stencil[1]])
    #fnum_i .= (fL .+ fR) ./ 2 - cache.cfl_loc./ 2 * (uR .- uL)
    #=
    @inbounds @simd for j in eachindex(fnum_i)
        fnum[i,j] = 0.5 * (fL[j] + fR[j]) - 0.5 * cache.cfl_loc * (uR[j] - uL[j])
    end
    =#
end