"""
     newton_lcp(M::AbstractMatrix, q::AbstractVector, z0::AbstractVector, w0::AbstractVector; maxiter::Int=1000, c=1)

Solve the linear complementarity problem LCP(`q`,`M`) with a Newton algorithm.
"""
function newton_lcp(M::AbstractMatrix, q::AbstractVector, z0::AbstractVector, w0::AbstractVector; maxiter::Int=1000, c=1)

    println("SOLVING WITH NEWTON LCP...")

    m = size(M)[1]
    Id = Matrix{eltype(M)}(I,m,m)

    # Initiate stopping criteria for the loop
    test = -1 # Striclty negative if In != Inprev, 0 else
    niter = 0

    # Initialize variables and indices
    z = copy(z0)
    w = copy(w0)
    In = w .- c .*z .<= 0
    Ac = .!In
    Inprev = zero(In)

    while niter < maxiter && test<0

        # First update the zero part
        z[Ac] .= 0
        w[In] .= 0

        # Solving a system for the remaining indices
        Mat = [M[:, In] -Id[:, Ac]]
        X = inv(Mat) * q

        z[In] .= X[1:sum(In)]
        w[Ac] .= X[sum(In) + 1 : end]

        Inprev .= In

        # Update indices
        @. In = w - c *z <= 0
        @. Ac = .!In

        # Testing stopping criterion
        if sum(In) == sum(Inprev)
            test = - norm((In .-Inprev) .> 0) + 1e-16
        else
            test = -1
        end

        niter += 1

    end

    if niter == maxiter
        println("Maximal number of iterations reached ("*string(maxiter)*").")
    else
        println(string(niter)*" iterations.")
    end

    z, w, niter

end