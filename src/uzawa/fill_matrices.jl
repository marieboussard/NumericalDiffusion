function fill_A!(A::AbstractMatrix, estimate::DiffEstimate)
    @unpack dt, uinit = estimate
    @unpack Nx, dx = estimate.params.mesh
    T = eltype(uinit)
    lamb = dt/dx
    for j in 1:Nx
        # A[mod1(j-1,Nx),j] = -one(T)*lamb
        # A[j,j] = one(T)*lamb
        # A[mod1(j-1,Nx)+Nx,j] = one(T)*lamb
        # A[j+Nx,j] = -one(T)*lamb
        A[mod1(j+1,Nx),j] = one(T)*lamb
        A[j,j] = -one(T)*lamb
        A[mod1(j+1,Nx)+Nx,j] = -one(T)*lamb
        A[j+Nx,j] = one(T)*lamb
    end
end

function fill_b!(b::AbstractVector, estimate::DiffEstimate)
    @unpack l, L = estimate
    @unpack Nx = estimate.params.mesh
    @views b_low = b[begin:Nx]
    @views b_up = b[Nx+1:end]
    b_low .= -l
    b_up .= L
    # for j in 1:Nx
    #     b_low[j] = -l[mod1(j+1,Nx)]
    #     b_up[j] = L[mod1(j+1,Nx)]
    # end
end

function fill_W!(W::AbstractMatrix, estimate::DiffEstimate, alpha::Real, beta::Real=1e-10)
    @unpack l, L = estimate
    @unpack Nx, dx = estimate.params.mesh
    for j in 1:Nx
        c = min(abs(L[j]-l[j])/dx, abs(L[mod1(j+1,Nx)] - l[mod1(j+1,Nx)])/dx)^alpha
        if c < beta
            W[j,j] = one(typeof(c))/beta#^alpha
        else
            W[j,j] = 1.0/c
            # W[j,j]=c
            # W[j,j] = one(typeof(c))
        end
    end
end