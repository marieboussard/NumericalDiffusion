struct UzawaSol{gtype<:AbstractVector, ptype<:AbstractVector}
    status::String
    gamma_opt::gtype
    p0::ptype
    popt::ptype
    mu::Float64
    niter::Int
    Jopt::Float64
    constraint_residual::Float64
    Gcgap::gtype  # |γopt - Gc|
    function UzawaSol(optimizer::Optimizer)
        @unpack W, gamma, Gc = optimizer
        Jopt = (W*(gamma .- Gc))'*(gamma .- Gc)
        Gcgap = gamma .- Gc
        # mul!(cache.Agamma, A, gamma)
        # constraint_residual = norm(max.(0.0,cache.Agamma .- b))
        if optimizer.iterate_gap <= optimizer.opts.eps && optimizer.constraint_residual <= optimizer.opts.eps_cons
            status = "SUCCESS"
            println("Convergence criteria reached! (eps<"*string(optimizer.opts.eps)*" and eps_cons<"*string(optimizer.opts.eps_cons)*") with "*string(optimizer.niter)*" iterations")
        elseif optimizer.niter == optimizer.opts.maxiter 
            status = "MAXITER"
            println("Stopped because maximum number of iterations was reached ("*string(optimizer.opts.maxiter)*")")
        else
            status = "FAILED"
            println("Stopped for unknow reason")
        end
        println("Constraint residual: "*string(optimizer.constraint_residual))
        new{typeof(optimizer.gamma), typeof(optimizer.p0)}(status, optimizer.gamma, optimizer.p0, optimizer.p, optimizer.mu, optimizer.niter, Jopt, optimizer.constraint_residual, Gcgap)
    end
end