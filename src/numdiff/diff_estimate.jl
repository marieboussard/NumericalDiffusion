struct DiffEstimate

    # PROBLEM COMPONENTS
    equation::Equation
    params::Parameters
    time_scheme::TimeScheme
    space_scheme::SpaceScheme

    # DATA
    u0
    u

    # QUANTIFICATION METHOD
    method::QuantifMethod
    qparams::qparamsType



    name::String
    function DiffEstimate(estimator::Estimator, name::String)
        new(name)
    end
end