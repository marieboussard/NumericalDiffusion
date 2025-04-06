function optimize_uzawa(Gc::AbstractVector, A::AbstractMatrix, b::AbstractVector; p0::AbstractVector=zero(b), gamma0::AbstractVector=zero(Gc), W::AbstractMatrix=Matrix{eltype(Gc)}(I,length(Gc),length(Gc)), maxiter::Int=10000, eps::Float64=1e-5, eps_cons::Float64=1e-5, start_with_gamma::Bool=false)
    println("SOLVING WITH UZAWA ALGORITHM...")
    optimizer = Optimizer(Gc, A, b; p0=p0, gamma0=gamma0, W=W, maxiter=maxiter, eps=eps, eps_cons=eps_cons, start_with_gamma=start_with_gamma)
    compute_mu!(optimizer)
    if start_with_gamma
        initstep!(optimizer)
    end
    while optimizer.niter < maxiter && (optimizer.iterate_gap > eps || optimizer.constraint_residual > eps_cons)
        performstep!(optimizer)
        loopfooter!(optimizer)
    end
    UzawaSol(optimizer)
end

function performstep!(optimizer::Optimizer)
    @unpack A, b, Gc, gamma, pprev, p, mu = optimizer
    @unpack Winv, Agamma, Atp, WAtp = optimizer.cache
    # 1 # Compute gamma as minimize the lagrangian without constraints
    mul!(Atp, A', pprev)
    mul!(WAtp, Winv, Atp)
    @. gamma = Gc - WAtp
    # 2 # Compute p as the orthogonal projection of the dual gradient descent iterate
    mul!(Agamma, A, gamma)
    for i in eachindex(p)
        p[i] = max(zero(eltype(pprev)), pprev[i]+mu*(Agamma[i]-b[i]))
    end
end

function compute_mu!(optimizer::Optimizer)
    @unpack W, A = optimizer
    optimizer.mu = minimum(diag(W))/opnorm(A)^2
end

function initstep!(optimizer::Optimizer)
    @unpack A, b, gamma0, gammaprev, p0, pprev, mu = optimizer
    @unpack Agamma = optimizer.cache
    mul!(Agamma, A, gamma0)
    for i in eachindex(p0)
        pprev[i] = max(zero(eltype(p0)), p0[i]+mu*(Agamma[i]-b[i]))
    end
    copyto!(gammaprev, gamma0)
end

function loopfooter!(optimizer::Optimizer)
    @unpack gammaprev, gamma, pprev, p, b = optimizer
    @unpack Agamma = optimizer.cache
    optimizer.iterate_gap = norm(gamma .- gammaprev)
    optimizer.constraint_residual = norm(max.(0.0, Agamma .- b))
    copyto!(gammaprev, gamma)
    copyto!(pprev, p)
    optimizer.niter += 1
end