
init_cache(::Euler) = EulerCache()
init_cache(::Rusanov) = RusanovCache()
# init_cache(::Rusanov, Nx::Int) = RusanovCache(Nx)