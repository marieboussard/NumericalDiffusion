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

function numflux!(::Rusanov, integrator::Integrator, u, i, args...)
    @unpack flux = integrator.equation
    @views uL = u[1,:]
    @views uR = u[2,:]
    @views intflux = integrator.flux[i,:]
    CFL_local!(integrator, u)
    intflux .= (flux(uL) .+ flux(uR)) ./ 2 - integrator.cache.cfl_loc./ 2 * (uR .- uL)
end