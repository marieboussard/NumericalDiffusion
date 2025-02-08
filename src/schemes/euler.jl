struct Euler <: TimeScheme end

struct EulerCache <: tcacheType end

compute_sL(::Euler, scheme::SpaceScheme) = compute_sL(scheme)
compute_sR(::Euler, scheme::SpaceScheme) = compute_sR(scheme)

numflux(::Euler, space_scheme::SpaceScheme, args...) = numflux(space_scheme, args...)

numflux!(::Euler, integrator::Integrator, args...) = numflux!(integrator.space_scheme, integrator, args...)