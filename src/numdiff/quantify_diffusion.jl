function quantify_diffusion(sol::Solution, method::QuantifMethod, params::QuantifParams, name=""; kwargs...)
    estimator = Estimator(data, method, params; kwargs...)
    perform_estimation!(estimator.method, estimator)
    DiffEstimate!(estimator, name)
end

function perform_estimation!(::Posteriori, estimator::Estimator)
end