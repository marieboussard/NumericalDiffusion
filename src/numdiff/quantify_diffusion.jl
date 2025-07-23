"""
    quantify_diffusion(sol::Solution, method::QuantifMethod, i=0; name="", <keyword arguments>)

Compute an estimation of the numerical diffusion produced by the numerical scheme used to obtain `sol`, using `method`.
"""
function quantify_diffusion(sol::Solution, method::QuantifMethod, i=0; name="", kwargs...)
    @btime estimator = Estimator(sol, method, i)
    @btime initialize_estimator!(estimator)
    @btime perform_estimation!(estimator.method, estimator; kwargs...)
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

"""
    compute_G_bounds!(estimator::Estimator)

Compute bounds on the numerical entropy flux G, that must be satisfied if G is a consistent entropic flux.
"""
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

