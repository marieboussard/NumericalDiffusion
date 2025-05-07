struct Euler <: TimeScheme end

struct EulerCache <: TimeCache end

get_sL(::Euler, scheme::SpaceScheme) = get_sL(scheme)
get_sR(::Euler, scheme::SpaceScheme) = get_sR(scheme)

# numflux(::Euler, space_scheme::SpaceScheme, args...) = numflux(space_scheme, args...)
# numflux!(::Euler, space_scheme::SpaceScheme, j::Int, args...) = numflux!(space_scheme, j, args...)

# numflux!(::Euler, integrator::Integrator, args...) = numflux!(integrator.space_scheme, integrator, args...)

### GLOBAL VERSION 

# Faster implementation with integrator
function numflux!(::Euler, integrator::Integrator)
    for j ∈ 1:integrator.params.mesh.Nx
        numflux!(integrator.space_scheme, integrator, j)
    end
    nothing
end

# General case
function numflux!(::Euler, space_scheme::SpaceScheme, params::Parameters, args...)
    for j ∈ 1:params.mesh.Nx
        numflux!(space_scheme, params, args..., j)
    end
    nothing
end

### LOCAL VERSION 
numflux!(::Euler, space_scheme::SpaceScheme, j::Int, i::Int, args...) = numflux!(space_scheme, j, i, args...)