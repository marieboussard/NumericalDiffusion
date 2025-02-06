init_cache(::Integrator) = IntegratorCache()

init_cache(::Euler) = EulerCache()
init_cache(::Rusanov) = RusanovCache()