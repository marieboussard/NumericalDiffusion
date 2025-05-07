## Caches depending on numerical schemes, used to compute uhat from utilde

# Euler scheme does not have any subcache, one needs to use the one in the estimator cache
linksubcache(::EulerCache, estimator::Estimator) = estimator.cache.subcache

# RK2 already has a subcache
linksubcache(cache::RK2Cache, ::Estimator) = cache.subcache


# Updating time cache

# EULER
function update_timecache!(::Euler, estimator::Estimator, args...)
    @unpack space_scheme, equation, cache = estimator
    update_subcache!(space_scheme, equation.dim, equation.eqtype, cache.subcache, equation, cache.utilde)
end

# RK2

function update_timecache!(::RK2, estimator::Estimator, time_cache::RK2Cache, u::AbstractVector)
    @unpack space_scheme, equation, cache = estimator
    @unpack u_bar, fcont_bar, fnum_bar, subcache = time_cache

    # 1 # update subcache for space scheme
    update_subcache!(space_scheme, equation.dim, equation.eqtype, subcache, equation, cache.utilde)

    # 2 # Apply space scheme with Euler to compute u_bar
    numflux!()

    # 3 # Compute fcont_bar = f(u_bar)
    flux!(equation.funcs, u_bar, fcont_bar)

end