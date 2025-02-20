# INIT CACHE FOR TIME SCHEMES
init_cache(::Euler) = EulerCache()

# INIT CACHE FOR SPACE SCHEMES
init_cache(::Rusanov) = RusanovCache()
init_cache(::Rusanov2D) = Rusanov2DCache()
init_cache(::Roe) = RoeCache()

# INIT CACHE FOR SOURCE
init_cache(ts::TopoSource, args...) = TopoSourceCache(ts, ts.source_discretize, args...)
init_cache(::NoSource, args...) = nothing