struct RusanovCache <: scacheType
    absDfcont
    function RusanovCache(Nx::Int)
        new(zeros(Float64, Nx))
    end
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

# function numflux!(::Rusanov, integrator::Integrator, u, i, args...)
#     @unpack flux = integrator.equation
#     @views uL = u[1,:]
#     @views uR = u[2,:]
#     @views intflux = integrator.flux[i,:]
#     CFL_local!(integrator, u)
#     intflux .= (flux(uL) .+ flux(uR)) ./ 2 - integrator.cache.cfl_loc./ 2 * (uR .- uL)
# end

function numflux!(::Rusanov, integrator::Integrator)
    @unpack equation, cache, params, uprev, fnum, fcont = integrator 
    @unpack cfl_loc = cache
    @unpack absDfcont = integrator.space_cache

    absDfcont .= abs.(equation.Dflux(uprev))
    
    @views uforward = uprev[2:end,:]
    @views ubackward = uprev[1:end-1,:]
    @views fforward = fcont[2:end,:]
    @views fbackward = fcont[1:end-1,:]
    # @views Dfforward = Dfcont[2:end,:]
    # @views Dfbackward = Dfcont[1:end-1,:]

    for i in 1:params.mesh.Nx-1
        cfl_loc[i] = max(absDfcont[i+1], absDfcont[i])
    end

    @. fnum[2:end-1,:] .= (fforward + fbackward)/2 - cfl_loc[1:end-1] /2 *(uforward - ubackward)

    # @show fnum[end,:]
    # @show (fcont[end,:] + fcont[1,:])/2 - max(absDfcont[end,:], absDfcont[1,:])/2*(uprev[end,:]-uprev[1,:])

    @. fnum[end,:] .= (fcont[end,:] + fcont[1,:])/2 - max(absDfcont[end,:], absDfcont[1,:])/2*(uprev[end,:]-uprev[1,:])

    #cache.cfl_loc = maximum(abs.(equation.Dflux(u)))

    # for i ∈ 2:integrator.params.mesh.Nx+1

    #     view_stencil!(integrator, i-1)
    #     # if equation.p==1
    #     #     uL = view(uprev,cache.stencil)[1,:]
    #     #     uR = view(uprev,cache.stencil)[2,:]
    #     # else
    #     #     uL = view(uprev,cache.stencil,:)[1,:]
    #     #     uR = view(uprev,cache.stencil,:)[2,:]
    #     # end
    #     CFL_local!(integrator, view(uprev,cache.stencil))
    #     fnum[i,:] .= (sum.(fcont[cache.stencil])) ./ 2 - cache.cfl_loc./ 2 * (view(uprev,cache.stencil)[2,:] .- view(uprev,cache.stencil)[1,:])
    #     #fnum[i,:] .= (flux(uL) .+ flux(uR)) ./ 2 - cache.cfl_loc./ 2 * (uR .- uL)
    # end
    fnum[1,:] .= fnum[end,:]
end