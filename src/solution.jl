struct Solution{utype<:AbstractArray}
    # PROBLEM COMPONENTS
    equation::Equation
    params::Parameters
    time_scheme::TimeScheme 
    space_scheme::SpaceScheme

    status::String      # SUCCESS or MAXITERS
    niter::Int          # Final number of iterations

    u::utype
    uinit::utype
    dt::Float64                  # Final timestep
    t::Float64                   # Time reached

    log::LogBook
    name::String

    function Solution(integrator::Integrator, name::String)
        if integrator.t == integrator.params.tf
            status = "SUCCES"
        elseif integrator.niter == integrator.opts.maxiter 
            status = "MAXITER"
        else
            status = "FAILED"
        end
        new{typeof(integrator.uinit)}(integrator.equation, integrator.params, integrator.time_scheme, integrator.space_scheme, status, integrator.niter, integrator.u, integrator.uinit, integrator.dt, integrator.t, integrator.log, name)
    end

end