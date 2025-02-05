struct Solution
    problem::Problem

    status::string      # SUCCESS or MAXITERS
    niter::Int          # Final number of iterations

    u
    dt                  # Final timestep
    t                   # Time reached

    log::LogBook
    name::string

end