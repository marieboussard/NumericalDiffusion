function utilde!(::SymmetricMD, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+1 : j+sL+2*sR, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    indices_short = view(indices, sL+sR+1:2*(sL+sR))
    @views ushorter = selectdim(uinit, 1, indices_short)
    computeK!(estimator.method.mdtype, ushorter, estimator.cache.K)
    @unpack K = estimator.cache
    if ndims(uinit)==1 # Scalar equations
        for i in 1:sL+sR
            utilde[i] = K[1]
        end
        for i in sL+sR+1:2*(sL+sR)
            utilde[i] = ushort[i]
        end
        for i in 2*(sL+sR)+1:3*(sL+sR)
            utilde[i] = K[1]
        end
    else    # Systems
        for r in 1:equation.p
            for i in 1:sL+sR
                utilde[i,r] = K[r]
            end
            for i in sL+sR+1:2*(sL+sR)
                utilde[i,r] = ushort[i,r]
            end
            for i in 2*(sL+sR)+1:3*(sL+sR)
                utilde[i,r] = K[r]
            end
        end
    end
end

function sourcetilde!(::SymmetricMD, estimator::Estimator, j::Int)# compute ̃zᵢʲ⁺¹/²
    @unpack znum = estimator.source_cache
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, sourceterm_tilde = estimator.cache
    # indices .= mod1.(j-2*sL-sR+1 : j+sL+2*sR, Nx)
    @views sourceshort = selectdim(znum, 1, indices)
    indices_short = view(indices, sL+sR+1:2*(sL+sR))
    @views sourceshorter = selectdim(znum, 1, indices_short)
    # computeK!(estimator.method.mdtype, sourceshorter, estimator.cache.S)
    S = computeK(estimator.method.mdtype, sourceshorter)
    @unpack S = estimator.cache
    for i in 1:sL+sR
        sourceterm_tilde[i] = S
    end
    for i in sL+sR+1:2*(sL+sR)
        sourceterm_tilde[i] = sourceshort[i]
    end
    for i in 2*(sL+sR)+1:3*(sL+sR)
        sourceterm_tilde[i] = S
    end
end

function init_bounds!(::SymmetricMD, estimator::Estimator, j::Int)
    @unpack equation, m, M, cache, entfun = estimator
    cache.GK = G(entfun, cache.K)
    if has_source(equation.source)
        cache.GK += Gsource(entfun, cache.K, cache.S)
    end
    m[j] = cache.GK
    M[j] = cache.GK
end

function update_bounds!(::SymmetricMD, estimator::Estimator, j::Int)
    @unpack equation, m, M, cache, entfun, dt = estimator
    @unpack sL, sR, utilde, uhat, eta_tilde, eta_hat, sourceterm_tilde = cache
    @unpack dx = estimator.params.mesh
    @views utilde_short = selectdim(utilde, 1, sL+1:3*sL+2*sR)
    eta!(entfun, utilde_short, eta_tilde)
    eta!(entfun, uhat, eta_hat)
    if has_source(equation.source)
        @views sourceterm_tilde_short = selectdim(sourceterm_tilde, 1, sL+1:3*sL+2*sR)
        etasource!(entfun, utilde_short, sourceterm_tilde_short, eta_tilde)
        etasource!(entfun, uhat, sourceterm_tilde_short, eta_hat)
    end
    for i=1:sL+sR
        M[j] += dx / dt *(eta_tilde[i] - eta_hat[i])
    end
    for i=sL+sR+1:2*(sL+sR)
        m[j] += dx / dt * (eta_hat[i] - eta_tilde[i])
    end
end