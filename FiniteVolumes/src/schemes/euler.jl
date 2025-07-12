struct Euler <: TimeScheme end

struct EulerCache <: TimeCache end

get_sL(::Euler, scheme::SpaceScheme) = get_sL(scheme)
get_sR(::Euler, scheme::SpaceScheme) = get_sR(scheme)

numflux(::Euler, space_scheme::SpaceScheme, args...) = numflux(space_scheme, args...)
numflux!(::Euler, space_scheme::SpaceScheme, j::Int, args...) = numflux!(space_scheme, j, args...)

numflux!(::Euler, integrator::Integrator, args...) = numflux!(integrator.space_scheme, integrator, args...)