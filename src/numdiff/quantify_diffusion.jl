function quantify_diffusion(sol::Solution, method::QuantifMethod; name="", kwargs...)
    estimator = Estimator(sol, method)
    initialize_estimator!(estimator)
    perform_estimation!(estimator.method, estimator; kwargs...)
    DiffEstimate(estimator, name)
end

function perform_estimation!(method::Posteriori, estimator::Estimator; kwargs...)
    @unpack Ginit, Gopt = estimator.method_cache
    compute_G_bounds!(estimator)
    @unpack m, M = estimator.method_cache
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
    @unpack mdtype, boundstype = estimator.method
    for j in 1:Nx
        utilde!(mdtype, boundstype, estimator, j)
        has_source(estimator.equation.source) ? sourcetilde!(mdtype, boundstype, estimator, j) : nothing
        uhat!(mdtype, boundstype, estimator)
        init_bounds!(mdtype, boundstype, estimator, j)
        update_bounds!(mdtype, boundstype, estimator, j)
    end
end

