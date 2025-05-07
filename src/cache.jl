# INIT CACHE FOR TIME SCHEMES
init_cache(::Euler, args...) = EulerCache()
init_cache(::RK2, args...) = RK2Cache(args...)

# INIT CACHE FOR SPACE SCHEMES
init_cache(::Rusanov, args...) = RusanovCache()
init_cache(::Rusanov2D, args...) = Rusanov2DCache()
init_cache(::Roe, args...) = RoeCache()
init_cache(scheme::HR, args...) = HRCache(scheme, args...)
init_cache(scheme::HR2D, args...) = HRCache(scheme, args...)

# INIT CACHE FOR SOURCE
init_cache(ts::TopoSource, args...) = TopoSourceCache(ts, ts.source_discretize, args...)
init_cache(::NoSource, args...) = nothing

# FILLING CACHES AFTER INTEGRATOR INITIALIZATION
fillcache!(::SpaceScheme, args...) = nothing

# SUBCACHE 
init_subcache(args...) = nothing