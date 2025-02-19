struct Euler <: TimeScheme end

struct EulerCache <: tcacheType end

get_sL(::Euler, scheme::SpaceScheme) = get_sL(scheme)
get_sR(::Euler, scheme::SpaceScheme) = get_sR(scheme)

numflux(::Euler, space_scheme::SpaceScheme, args...) = numflux(space_scheme, args...)

numflux!(::Euler, integrator::Integrator, args...) = numflux!(integrator.space_scheme, integrator, args...)