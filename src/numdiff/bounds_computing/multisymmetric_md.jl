init_indices(::SymmetricMD, ::MultiBounds, sL::Int, sR::Int) = zeros(Int64, 3*(sL+sR)+1)

# ONE DIMENSIONAL SCALAR EQUATIONS
init_utilde(::SymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR)+1)
init_uhat(::SymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)+1)
init_ftilde(::SymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)+2)

function utilde!(::SymmetricMD, ::MultiBounds, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR : j+sL+2*sR, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    indices_short = view(indices, sL+sR+1:2*(sL+sR)+1)
    @views ushorter = selectdim(uinit, 1, indices_short)
    computeK!(estimator.method.mdtype, ushorter, estimator.cache.mdcache.K)
    @unpack K = estimator.cache.mdcache
    if ndims(uinit)==1 # Scalar equations
        for i in 1:sL+sR
            utilde[i] = K[1]
        end
        for i in sL+sR+1:2*(sL+sR)+1
            utilde[i] = ushort[i]
        end
        for i in 2*(sL+sR)+2:3*(sL+sR)+1
            utilde[i] = K[1]
        end
    # else    # Systems
    #     for r in 1:equation.p
    #         for i in 1:sL+sR
    #             utilde[i,r] = K[r]
    #         end
    #         for i in sL+sR+1:2*(sL+sR)
    #             utilde[i,r] = ushort[i,r]
    #         end
    #         for i in 2*(sL+sR)+1:3*(sL+sR)
    #             utilde[i,r] = K[r]
    #         end
    #     end
    end
end

function uhat!(::SymmetricMD, ::MultiBounds, estimator::Estimator)
    @unpack time_scheme, space_scheme, params, equation, cache, space_cache, source_cache, dt = estimator
    @unpack sL, sR, utilde, ftilde, fcont_tilde, uhat = cache
    @unpack dx = params.mesh
    flux!(equation.funcs, utilde, fcont_tilde)
    update_cflcache!(equation.dim, equation.eqtype, equation.funcs, utilde, cache.cfl_cache)
    for i in 1:2*(sL+sR)+2
        numflux!(time_scheme, space_scheme, i+sL-1, params, equation, cache, space_cache, ftilde, fcont_tilde, utilde, i)
    end
    for i in 1:2*(sL+sR)+1
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

# function init_bounds!(::SymmetricMD, ::DefaultBounds, estimator::Estimator, j::Int)
#     @unpack m, M = estimator
#     # cache.mdcache.GK = G(entfun, cache.mdcache.K)
#     # if has_source(equation.source)
#     #     cache.mdcache.GK += Gsource(entfun, cache.mdcache.K, cache.mdcache.S)
#     # end
#     # m[j] = cache.mdcache.GK
#     # M[j] = cache.mdcache.GK
#     m[j] = zero(eltype(m))
#     M[j] = zero(eltype(M))
# end