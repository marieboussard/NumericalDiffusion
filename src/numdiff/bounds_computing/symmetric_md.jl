mutable struct SymmetricMDCache{ktype<:AbstractArray} <: ModifiedDataCache
    K::ktype
    S::Float64
    GK::Float64
    function SymmetricMDCache(equation::Equation, u::AbstractArray)
        K = zeros(eltype(u), equation.p)
        new{typeof(K)}(K, zero(Float64), zero(Float64))
    end
end

init_cache(::SymmetricMD, equation::Equation, u::AbstractArray) = SymmetricMDCache(equation, u)

# ONE DIMENSIONAL SCALAR EQUATIONS
init_utilde(::SymmetricMD, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR))
init_uhat(::SymmetricMD, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR))
init_ftilde(::SymmetricMD, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)+1)
# init_K(::OneD, ::Scalar) = zero(Float64)

# ONE DIMENSIONAL SYSTEMS
init_utilde(::SymmetricMD, ::OneD, ::System, u::Matrix{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR), size(u)[2])
init_uhat(::SymmetricMD, ::OneD, ::System, u::Matrix{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR), size(u)[2])
init_ftilde(::SymmetricMD, ::OneD, ::System, u::Matrix{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)+1, size(u)[2])
# init_K(::OneD, ::Scalar) = zeros(Float64, size(u)[2])


function utilde!(::SymmetricMD, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+1 : j+sL+2*sR, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    indices_short = view(indices, sL+sR+1:2*(sL+sR))
    @views ushorter = selectdim(uinit, 1, indices_short)
    computeK!(estimator.method.mdtype, ushorter, estimator.cache.mdcache.K)
    @unpack K = estimator.cache.mdcache
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
    estimator.cache.mdcache.S = computeK(estimator.method.mdtype, sourceshorter)
    @unpack S = estimator.cache.mdcache
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

function uhat!(::SymmetricMD, estimator::Estimator)
    @unpack time_scheme, space_scheme, params, equation, cache, space_cache, source_cache, dt = estimator
    @unpack sL, sR, utilde, ftilde, fcont_tilde, uhat = cache
    @unpack dx = params.mesh
    flux!(equation.funcs, utilde, fcont_tilde)
    update_cflcache!(equation.dim, equation.eqtype, equation.funcs, utilde, cache.cfl_cache)
    for i in 1:2*(sL+sR)+1
        numflux!(time_scheme, space_scheme, i+sL-1, params, equation, cache, space_cache, ftilde, fcont_tilde, utilde, i)
    end
    for i in 1:2*(sL+sR)
        for r in 1:equation.p
            uhat[i,r] = utilde[i+sL,r] - dt/dx * (ftilde[i+1,r] - ftilde[i,r])
        end
    end
    if has_source(equation.source)
        discretize_sourceterm!(equation.dim, equation.source.source_discretize, utilde, cache.sourceterm_tilde, source_cache)
        for i in 1:2*(sL+sR)
            for j in 1:equation.p
                uhat[i,j] += dt * cache.sourceterm_tilde[i,j]
            end
        end
    end
end

function init_bounds!(::SymmetricMD, estimator::Estimator, j::Int)
    @unpack equation, m, M, cache, entfun = estimator
    cache.mdcache.GK = G(entfun, cache.mdcache.K)
    if has_source(equation.source)
        cache.mdcache.GK += Gsource(entfun, cache.mdcache.K, cache.mdcache.S)
    end
    m[j] = cache.mdcache.GK
    M[j] = cache.mdcache.GK
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