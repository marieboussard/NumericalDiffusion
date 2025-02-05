struct Solution
    problem::Problem
    parameters::Parameters

    status::string      # SUCCESS or IN PROGRESS or MAXITERS
    niter::Int          # Final number of iterations

    u
    dt                  # Final timestep
    t                   # Time reached

    log::LogBook
    name::string

end