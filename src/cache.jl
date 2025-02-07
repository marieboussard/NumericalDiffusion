init_cache(sL, sR) = IntegratorCache(sL, sR)

init_cache(::Euler) = EulerCache()
init_cache(::Rusanov) = RusanovCache()