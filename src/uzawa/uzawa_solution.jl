struct UzawaSol{gtype<:AbstractVector, ptype<:AbstractVector}
    status::String
    gamma_opt::gtype
    p0::ptype
    popt::ptype
    mu::Float64
    niter::Int
    constraint_residual::Float64
    function UzawaSol(optimizer::Optimizer)
        @unpack A, b, gamma, cache = optimizer
        mul!(cache.Agamma, A, gamma)
        constraint_residual = norm(cache.Agamma .- b)
        if optimizer.iterate_gap <= optimizer.opts.eps
            status = "SUCCESS"
            println("Convergence criteria reached! (eps<"*string(optimizer.opts.eps)*")")
        elseif optimizer.niter == optimizer.opts.maxiter 
            status = "MAXITER"
            println("Stopped because maximum number of iterations was reached ("*string(optimizer.opts.maxiter)*")")
        else
            status = "FAILED"
            println("Stopped for unknow reason")
        end
        println("Constraint residual: "*string(constraint_residual))
        new{typeof(optimizer.gamma), typeof(optimizer.p0)}(status, optimizer.gamma, optimizer.p0, optimizer.p, optimizer.mu, optimizer.niter, constraint_residual)
    end
end