function quantify_diffusion(sol::Solution, method::QuantifMethod, name=""; kwargs...)
    estimator = Estimator(sol, method)
    perform_estimation!(estimator.method, estimator; kwargs...)
    DiffEstimate(estimator, name)
end

function perform_estimation!(method::Posteriori, estimator::Estimator; kwargs...)
    @unpack Ginit, Gopt = estimator.method_cache

    eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
    eta!(estimator.entfun, estimator.u, estimator.etacont)

    compute_G_bounds!(estimator)
    @unpack m, M = estimator
    for i in eachindex(Ginit)
        Ginit[i] = 0.5*(m[i]+M[i])
    end
    optsol = optimize(gamma -> J(estimator, gamma), Ginit; kwargs...)
    copyto!(Gopt, Optim.minimizer(optsol))
    estimator.method_cache.Jopt = Optim.minimum(optsol)
    diffusion!(method, estimator)
end

function compute_G_bounds!(estimator::Estimator)
    @unpack Nx = estimator.params.mesh
    @unpack mdtype = estimator.method
    # @show estimator.uinit
    for j in 1:Nx
        utilde!(mdtype, estimator, j)
        # @show estimator.cache.utilde
        uhat!(estimator)
        # @show estimator.cache.uhat
        init_bounds!(mdtype, estimator, j)
        # @show estimator.m, estimator.M
        update_bounds!(mdtype, estimator, j)
        # @show estimator.m, estimator.M
    end
end

function uhat!(estimator::Estimator)
    @unpack time_scheme, space_scheme, params, equation, cache, space_cache, dt = estimator
    @unpack sL, sR, utilde, ftilde, fcont_tilde, uhat = cache
    @unpack dx = params.mesh
    flux!(equation.funcs, utilde, fcont_tilde)
    for i in 1:2*(sL+sR)+1
        numflux!(time_scheme, space_scheme, i+sL-1, params, equation, cache, space_cache, ftilde, fcont_tilde, utilde, i)
    end
    for i in 1:2*(sL+sR)
        # uhat[i] = utilde[i+sL] - dt/dx * (ftilde[i] - ftilde[mod1(i-1, 2*(sL+sR))])
        uhat[i] = utilde[i+sL] - dt/dx * (ftilde[i+1] - ftilde[i])
    end
end