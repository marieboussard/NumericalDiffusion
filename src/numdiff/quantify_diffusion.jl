function quantify_diffusion(sol::Solution, method::QuantifMethod, params::QuantifParams, name=""; kwargs...)
    estimator = Estimator(data, method, params)
    perform_estimation!(estimator.method, estimator; kwargs...)
    DiffEstimate!(estimator, name)
end

function perform_estimation!(::Posteriori, estimator::Estimator; kwargs...)
    compute_G_bounds!(estimator)
    optsol = optimize(gamma -> J(estimator, gamma), Ginit; kwargs...)
    estimator.Gopt, estimator.Jopt = Optim.minimizer(optsol), Optim.minimum(optsol)
    diffusion!(estimator)
end