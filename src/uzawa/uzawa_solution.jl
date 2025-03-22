struct UzawaSol{gtype<:AbstractVector, ptype<:AbstractVector}
    status::String
    gamma_opt::gtype
    p0::ptype
    popt::ptype
    mu::Float64
    niter::Int
    function UzawaSol(optimizer::Optimizer)
        if optimizer.iterate_gap <= optimizer.opts.eps
            status = "SUCCES"
        elseif integrator.niter == integrator.opts.maxiter 
            status = "MAXITER"
        else
            status = "FAILED"
        end
        new{typeof(optimizer.gamma), typeof(optimizer.p0)}(status, optimizer.gamma, optimizer.p0, optimizer.p, optimizer.mu, optimizer.niter)
    end
end