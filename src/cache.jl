init_cache() = IntegratorCache()

init_cache(::Euler) = EulerCache()
init_cache(::Rusanov) = RusanovCache()