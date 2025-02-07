struct RusanovCache <: Cache
    sL::Int
    sR::Int
    
    function RusanovCache()
        new(1, 1)
    end
end

struct Rusanov <: SpaceScheme
    cache::RusanovCache
    Rusanov()=new(RusanovCache())
end

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
    @unpack equation, cache, params, uprev = integrator
    @unpack flux = equation

    for i ∈ 2:integrator.params.mesh.Nx+1

        view_stencil!(integrator, i-1)
        if equation.p==1
            uL = view(uprev,cache.stencil)[1,:]
            uR = view(uprev,cache.stencil)[2,:]
        else
            uL = view(uprev,cache.stencil,:)[1,:]
            uR = view(uprev,cache.stencil,:)[2,:]
        end
        CFL_local!(integrator, [uL; uR])
        integrator.flux[i,:] .= (flux(uL) .+ flux(uR)) ./ 2 - cache.cfl_loc./ 2 * (uR .- uL)
    end
    integrator.flux[1,:] .= integrator.flux[end,:]
end