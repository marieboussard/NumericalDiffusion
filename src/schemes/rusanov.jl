struct Rusanov <: SpaceScheme
    cache::RusanovCache
    Rusanov()=new(RusanovCache())
end

struct RusanovCache <: Cache
    sL::Int
    sR::Int
    A
    
    function RusanovCache()
        new(1, 1, 0.0)
    end
end

compute_sL(::Rusanov) = 1
compute_sR(::Rusanov) = 1

function numflux(scheme::Rusanov, integrator::Integrator, u, args...)
    uL, uR = u[1,:], u[2,:]
    scheme.cache.A = CFL_cond(integrator, u)
    (flux(equation, uL) .+ flux(equation, uR)) / 2 .- scheme.cache.A / 2 * (uR .- uL)
end