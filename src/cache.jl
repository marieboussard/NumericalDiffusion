# INIT CACHE FOR SPACE SCHEMES
init_cache(::Rusanov, u::AbstractArray, dim::EquationDim, args...) = RusanovCache(u)
init_cache(::Rusanov2D, args...) = Rusanov2DCache()
init_cache(::Roe, args...) = RoeCache()
init_cache(scheme::HR, args...) = HRCache(scheme, args...)
init_cache(scheme::HR2D, args...) = HRCache(scheme, args...)

# INIT CACHE FOR TIME SCHEMES
# init_cache(::Euler, space_cache::SpaceCache, args...) = EulerCache(space_cache, args...)


function init_cache(::Euler, space_scheme::SpaceScheme, args...)
    space_cache = init_cache(space_scheme, args...)
    EulerCache(space_cache)
end

init_cache(::RK2, space_scheme, u::AbstractArray, dim::EquationDim, params::Parameters, dt::Float64, args...) = RK2Cache(space_scheme, u, dim, params, dt)

# INIT CACHE FOR SOURCE
init_cache(ts::TopoSource, args...) = TopoSourceCache(ts, ts.source_discretize, args...)
init_cache(::NoSource, args...) = nothing

# FILLING CACHES AFTER INTEGRATOR INITIALIZATION
fillcache!(::SpaceScheme, args...) = nothing