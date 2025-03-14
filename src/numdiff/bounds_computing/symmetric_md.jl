function utilde!(::SymmetricMD, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+1 : j+sL+2*sR, Nx)
    @views ushort = uinit[indices]
    indices_short = view(indices, sL+sR+1:2*(sL+sR))
    @views ushorter = uinit[indices_short]
    # @show @allocated ushorter = view(ushort, sL+sR+1:2*(sL+sR))
    # @show @allocated ushorter = ushort[sL+sR+1:2*(sL+sR)]
    estimator.cache.K = computeK(estimator.method.mdtype, ushorter)
    @unpack K = estimator.cache
    for i in 1:sL+sR
        utilde[i] = K
    end
    for i in sL+sR+1:2*(sL+sR)
        utilde[i] = ushort[i]
    end
    for i in 2*(sL+sR)+1:3*(sL+sR)
        utilde[i] = K
    end
end

function init_bounds!(::SymmetricMD, estimator::Estimator, j::Int)
    @unpack m, M, cache, entfun = estimator
    cache.GK = G(entfun, cache.K)
    m[j] = cache.GK
    M[j] = cache.GK
end

function update_bounds!(::SymmetricMD, estimator::Estimator, j::Int)
    @unpack m, M, cache, entfun, dt = estimator
    @unpack sL, sR, utilde, uhat, eta_tilde, eta_hat = cache
    @unpack dx = estimator.params.mesh
    eta!(entfun, view(utilde, sL+1:3*sL+2*sR), eta_tilde)
    eta!(entfun, uhat, eta_hat)
    for i=1:sL+sR
        M[j] += dx / dt *(eta_tilde[i] - eta_hat[i])
    end
    for i=sL+sR+1:2*(sL+sR)
        m[j] += dx / dt * (eta_hat[i] - eta_tilde[i])
    end
end