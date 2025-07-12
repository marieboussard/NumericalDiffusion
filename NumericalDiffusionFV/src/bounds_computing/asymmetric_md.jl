mutable struct AsymmetricMDCache <: ModifiedDataCache end

init_cache(::AsymmetricMD, equation::Equation, u::AbstractArray) = AsymmetricMDCache()
init_indices(::AsymmetricMD, ::DefaultBounds, sL::Int, sR::Int) = zeros(Int64, 3*(sL+sR)-2)

# ONE DIMENSIONAL SCALAR EQUATIONS
init_utilde(::AsymmetricMD, ::DefaultBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR)-2)
init_uhat(::AsymmetricMD, ::DefaultBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)-2)
init_ftilde(::AsymmetricMD, ::DefaultBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)-1)


function utilde!(::AsymmetricMD, ::DefaultBounds, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+2 : j+sL+2*sR-1, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    if ndims(uinit)==1 # Scalar equations
        for i in 1:sL+sR-1
            utilde[i] = ushort[sL+sR]
        end
        for i in sL+sR:2*(sL+sR)-1
            utilde[i] = ushort[i]
        end
        for i in 2*(sL+sR):3*(sL+sR)-2
            utilde[i] = ushort[2*(sL+sR)-1]
        end
    # else    # Systems
    #     for r in 1:equation.p
    #         for i in 1:sL+sR-1
    #             utilde[i,r] = K[r]
    #         end
    #         for i in sL+sR:2*(sL+sR)-1
    #             utilde[i,r] = ushort[i,r]
    #         end
    #         for i in 2*(sL+sR):3*(sL+sR)-2
    #             utilde[i,r] = K[r]
    #         end
    #     end
    end
end

function uhat!(::AsymmetricMD, ::DefaultBounds, estimator::Estimator)
    @unpack time_scheme, space_scheme, params, equation, cache, space_cache, source_cache, dt = estimator
    @unpack sL, sR, utilde, ftilde, fcont_tilde, uhat = cache
    @unpack dx = params.mesh
    flux!(equation.funcs, utilde, fcont_tilde)
    update_cflcache!(equation.dim, equation.eqtype, equation.funcs, utilde, cache.cfl_cache)
    for i in 1:2*(sL+sR)-1
        numflux!(time_scheme, space_scheme, i+sL-1, params, equation, cache, space_cache, ftilde, fcont_tilde, utilde, i)
    end
    for i in 1:2*(sL+sR)-2
        for r in 1:equation.p
            uhat[i,r] = utilde[i+sL,r] - dt/dx * (ftilde[i+1,r] - ftilde[i,r])
        end
    end
    # if has_source(equation.source)
    #     discretize_sourceterm!(equation.dim, equation.source.source_discretize, utilde, cache.sourceterm_tilde, source_cache)
    #     for i in 1:2*(sL+sR)
    #         for j in 1:equation.p
    #             uhat[i,j] += dt * cache.sourceterm_tilde[i,j]
    #         end
    #     end
    # end
end

function init_bounds!(::AsymmetricMD, ::DefaultBounds, estimator::Estimator, j::Int)
    @unpack equation, cache, entfun, uinit = estimator
    @unpack m, M = estimator.method_cache
    @unpack sL, sR = cache
    @unpack Nx = estimator.params.mesh
    # if has_source(equation.source)
    #     cache.mdcache.GK += Gsource(entfun, cache.mdcache.K, cache.mdcache.S)
    # end
    M[j] = G(entfun, uinit[mod1(j-sL+1, Nx)])
    m[j] = G(entfun, uinit[mod1(j+sR, Nx)])
end

function update_bounds!(::AsymmetricMD, ::DefaultBounds, estimator::Estimator, j::Int)
    @unpack equation, cache, entfun, dt = estimator
    @unpack m, M = estimator.method_cache
    @unpack sL, sR, utilde, uhat, eta_tilde, eta_hat, sourceterm_tilde = cache
    @unpack dx = estimator.params.mesh
    @views utilde_short = selectdim(utilde, 1, sL+1:3*sL+2*sR-2)
    eta!(entfun, utilde_short, eta_tilde)
    eta!(entfun, uhat, eta_hat)
    # if has_source(equation.source)
    #     @views sourceterm_tilde_short = selectdim(sourceterm_tilde, 1, sL+1:3*sL+2*sR)
    #     etasource!(entfun, utilde_short, sourceterm_tilde_short, eta_tilde)
    #     etasource!(entfun, uhat, sourceterm_tilde_short, eta_hat)
    # end
    for i=1:sL+sR-1
        M[j] += dx / dt *(eta_tilde[i] - eta_hat[i])
    end
    for i=sL+sR:2*(sL+sR)-2
        m[j] += dx / dt * (eta_hat[i] - eta_tilde[i])
    end
end