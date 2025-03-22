function optimize_uzawa(Gc::AbstractVector, A::AbstractMatrix, b::AbstractVector; p0::AbstractVector=zero(b), W::AbstractMatrix=Matrix{eltype(Gc)}(I,length(Gc),length(Gc)), maxiter::Int=100, eps::Float64=1e-3)
    optimizer = Optimizer(Gc, A, b; p0=p0, W=W, maxiter=maxiter, eps=eps)
    compute_mu!(optimizer)
    while optimizer.niter < maxiter && optimizer.iterate_gap > eps
        performstep!(optimizer)
        loopfooter!(optimizer)
    end
    UzawaSol(optimizer)
end

function performstep!(optimizer::Optimizer)
    @unpack W, A, b, Gc, gamma, pprev, p, mu = optimizer
    @unpack Agamma, Atp, WAtp = optimizer.cache
    # 1 # Compute gamma as minimize the lagrangian without constraints
    mul!(Atp, A', pprev)
    mul!(WAtp, W, Atp)
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

function loopfooter!(optimizer::Optimizer)
    @unpack gammaprev, gamma, pprev, p = optimizer
    optimizer.iterate_gap = norm(gamma .- gammaprev)
    copyto!(gammaprev, gamma)
    copyto!(pprev, p)
    optimizer.niter += 1
end