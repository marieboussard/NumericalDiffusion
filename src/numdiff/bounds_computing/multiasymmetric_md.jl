init_indices(::AsymmetricMD, ::MultiBounds, sL::Int, sR::Int) = zeros(Int64, 3*(sL+sR)-1)

# ONE DIMENSIONAL SCALAR EQUATIONS
init_utilde(::AsymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 3*(sL+sR)-1)
init_uhat(::AsymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR)-1)
init_ftilde(::AsymmetricMD, ::MultiBounds, ::OneD, ::Scalar, u::Vector{Float64}, sL::Int, sR::Int) = zeros(eltype(u), 2*(sL+sR))

function utilde!(::AsymmetricMD, ::MultiBounds, estimator::Estimator, j::Int)# compute ̃uᵢʲ⁺¹/²
    @unpack uinit = estimator
    @unpack Nx = estimator.params.mesh
    @unpack sL, sR, indices, utilde = estimator.cache
    indices .= mod1.(j-2*sL-sR+1 : j+sL+2*sR-1, Nx)
    @views ushort = selectdim(uinit, 1, indices)
    if ndims(uinit)==1 # Scalar equations
        for i in 1:sL+sR-1
            utilde[i] = ushort[sL+sR]
        end
        for i in sL+sR:2*(sL+sR)
            utilde[i] = ushort[i]
        end
        for i in 2*(sL+sR)+1:3*(sL+sR)-1
            utilde[i] = ushort[2*(sL+sR)]
        end
    end
end

function uhat!(::AsymmetricMD, ::MultiBounds, estimator::Estimator)
    # @unpack time_scheme, space_scheme, params, equation, cache, space_cache, source_cache, dt = estimator
    # @unpack sL, sR, utilde, ftilde, fcont_tilde, uhat = cache
    # @unpack dx = params.mesh
    # flux!(equation.funcs, utilde, fcont_tilde)
    # update_cflcache!(equation.dim, equation.eqtype, equation.funcs, utilde, cache.cfl_cache)
    # for i in 1:2*(sL+sR)
    #     numflux!(time_scheme, space_scheme, i+sL-1, params, equation, cache, space_cache, ftilde, fcont_tilde, utilde, i)
    # end
    # for i in 1:2*(sL+sR)-1
    #     for r in 1:equation.p
    #         uhat[i,r] = utilde[i+sL,r] - dt/dx * (ftilde[i+1,r] - ftilde[i,r])
    #     end
    # end

    @unpack time_scheme, space_scheme, time_cache, params, equation, cache, source_cache, dt = estimator
    @unpack sL, sR, utilde, ftilde, uhat = cache
    @unpack dx = params.mesh

    global_numflux!(time_scheme, space_scheme, time_cache, equation, utilde, ftilde, sL, 3*sL+2*sR-1, params.mesh.Nx, -sL+1)
    for i in 1:2*(sL+sR)-1
        for r in 1:equation.p
            uhat[i,r] = utilde[i+sL,r] - dt/dx * (ftilde[i+1,r] - ftilde[i,r])
        end
    end

end

# BOUNDS

function compute_DG_bounds!(::AsymmetricMD, estimator::Estimator)
    @unpack entfun, uinit, etacont, etacont_init, dt = estimator
    @unpack Nx, dx = estimator.params.mesh
    @unpack mdtype, boundstype = estimator.method
    @unpack l, L = estimator.method_cache
    @unpack sL, sR, eta_tilde, eta_hat, utilde, uhat = estimator.cache
    for j in 1:Nx
        L[j] = etacont_init[j] - etacont[j]
        l[j] = dt / dx * (G(entfun, uinit[mod1(j+sR, Nx)]) - G(entfun, uinit[mod1(j-sL, Nx)]))
        utilde!(mdtype, boundstype, estimator, j)
        uhat!(mdtype, boundstype, estimator)
        @views utilde_short = selectdim(utilde, 1, sL+1:3*sL+2*sR-1)
        eta!(entfun, utilde_short, eta_tilde)
        eta!(entfun, uhat, eta_hat)
        for i in 1:sL+sR-1
            l[j] += eta_hat[i] - eta_tilde[i]
        end
        for i in sL+sR+1:2*(sL+sR)-1
            l[j] += eta_hat[i] - eta_tilde[i]
        end
    end
end