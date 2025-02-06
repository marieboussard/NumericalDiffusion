struct Euler <: TimeScheme end

struct EulerCache <: Cache
end

compute_sL(::Euler, scheme::space_scheme) = compute_sL(scheme)
compute_sR(::Euler, scheme::space_scheme) = compute_sR(scheme)


numflux(::Euler, integrator::Integrator, args...) = numflux(integrator.space_scheme, integrator, args...)

num_flux