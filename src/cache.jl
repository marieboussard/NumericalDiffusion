# INIT CACHE FOR TIME SCHEMES
init_cache(::Euler) = EulerCache()

# INIT CACHE FOR SPACE SCHEMES
init_cache(::Rusanov) = RusanovCache()

# INIT CACHE FOR SOURCE
init_cache(ts::TopoSource, args...) = TopoSourceCache(ts, args...)
init_cache(::NoSource, args...) = nothing