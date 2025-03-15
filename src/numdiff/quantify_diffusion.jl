function quantify_diffusion(sol::Solution, method::QuantifMethod, name=""; kwargs...)
    estimator = Estimator(sol, method)
    perform_estimation!(estimator.method, estimator; kwargs...)
    DiffEstimate(estimator, name)
end

function perform_estimation!(method::Posteriori, estimator::Estimator; kwargs...)
    @unpack Ginit, Gopt = estimator.method_cache
    eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
    eta!(estimator.entfun, estimator.u, estimator.etacont)
    if has_source(estimator.equation.source)
        etasource!(estimator.entfun, estimator.uinit, estimator.source_cache, estimator.etacont_init)
        etasource!(estimator.entfun, estimator.u, estimator.source_cache, estimator.etacont)
    end

    compute_G_bounds!(estimator)
    @unpack m, M = estimator
    for i in eachindex(Ginit)
        Ginit[i] = 0.5*(m[i]+M[i])
    end
    optsol = optimize(gamma -> J(estimator, gamma), Ginit; g_tol=1e-20, iterations=200000, method=LBFGS(), autodiff=:forward, kwargs...)
    copyto!(Gopt, Optim.minimizer(optsol))
    estimator.method_cache.Jopt = Optim.minimum(optsol)

    diffusion!(method, estimator)
end

function compute_G_bounds!(estimator::Estimator)
    @unpack Nx = estimator.params.mesh
    @unpack mdtype = estimator.method
    for j in 1:Nx
        utilde!(mdtype, estimator, j)
        has_source(estimator.equation.source) ? sourcetilde!(mdtype, estimator, j) : nothing
        uhat!(estimator)
        init_bounds!(mdtype, estimator, j)
        update_bounds!(mdtype, estimator, j)
    end
end

function uhat!(estimator::Estimator)
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