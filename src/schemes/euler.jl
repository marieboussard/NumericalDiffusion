struct Euler <: TimeScheme end

struct EulerCache{sctype<:SpaceCache} <: TimeCache 
    space_cache::sctype
end

get_sL(::Euler, scheme::SpaceScheme) = get_sL(scheme)
get_sR(::Euler, scheme::SpaceScheme) = get_sR(scheme)

get_name(::Euler) = "Euler" 

# numflux(::Euler, space_scheme::SpaceScheme, args...) = numflux(space_scheme, args...)
# numflux!(::Euler, space_scheme::SpaceScheme, j::Int, args...) = numflux!(space_scheme, j, args...)

# """
#     numflux!(::Euler, integrator::Integrator, args...)
# """
# numflux!(::Euler, integrator::Integrator, args...) = numflux!(integrator.space_scheme, integrator, args...)


function global_numflux!(::Euler, space_scheme::SpaceScheme, time_cache::TimeCache, equation::Equation, u::AbstractArray, fnum::AbstractArray, jstart::Int=1, jend::Int=size(u)[1], Nx=size(u)[1], shift::Int=0)
    # Global computation of necessary quantities
    update_cache!(time_cache.space_cache, u, equation.funcs)#, jstart, jend)
    # Iterate on each interface to compute the flux
    for j ∈ jstart:jend
        numflux!(space_scheme, time_cache.space_cache, equation, u, fnum, j, Nx, j+shift)
    end
    nothing
end

global_numflux!(ts::Euler, integrator::Integrator) = global_numflux!(ts, integrator.space_scheme, integrator.time_cache, integrator.equation, integrator.uprev, integrator.fnum)