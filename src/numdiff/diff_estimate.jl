struct DiffEstimate
    name::String
    function DiffEstimate(estimator::Estimator, name::String)
        new(name)
    end
end