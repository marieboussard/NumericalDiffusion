abstract type AbstractNormWeights end
struct AlphaWeights{T<:Real} <: AbstractNormWeights
    alpha::T
end
struct MinimizedAlphaWeights{T<:Real} <: AbstractNormWeights
    alpha::T
end
struct AbsWeights{T<:Real} <: AbstractNormWeights
    alpha::T
end
AbsWeights()=new{Int}(1)

function fill_A!(A::AbstractMatrix, estimate::DiffEstimate)
    @unpack dt, uinit = estimate
    @unpack Nx, dx = estimate.params.mesh
    T = eltype(uinit)
    lamb = dt/dx
    for j in 1:Nx
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
end

function fill_W!(aw::AbsWeights, W::AbstractMatrix, estimate::DiffEstimate, beta::Real=1e-10)
    @unpack alpha = aw
    @unpack l, L = estimate
    @unpack Nx, dx = estimate.params.mesh
    @unpack uinit = estimate
    for j in 1:Nx
        c = abs(uinit[j] - uinit[mod1(j+1,Nx)])^alpha
        if c < beta
            W[j,j] = one(typeof(c))/beta
        else
            W[j,j] = 1.0/c
        end
    end
end

function fill_W!(aw::AlphaWeights, W::AbstractMatrix, estimate::DiffEstimate, beta::Real=1e-10)
    @unpack alpha = aw
    @unpack l, L = estimate
    @unpack Nx, dx = estimate.params.mesh
    @unpack uinit = estimate
    for j in 1:Nx
        c = min(abs(L[j]-l[j])/dx, abs(L[mod1(j+1,Nx)] - l[mod1(j+1,Nx)])/dx)^alpha
        if c < beta
            W[j,j] = one(typeof(c))/beta
        else
            W[j,j] = 1.0/c
        end
    end
end

function fill_W!(aw::MinimizedAlphaWeights, W::AbstractMatrix, estimate::DiffEstimate, beta::Real=1e-10, betamax::Real = 5)
    @unpack alpha = aw
    @unpack l, L = estimate
    @unpack Nx, dx = estimate.params.mesh
    @unpack uinit = estimate
    for j in 1:Nx
        c = min(abs(L[j]-l[j])/dx, abs(L[mod1(j+1,Nx)] - l[mod1(j+1,Nx)])/dx)^alpha
        if c < beta
            W[j,j] = one(typeof(c))/beta
        elseif c > betamax
            println("max treshold reached")
            W[j,j] = one(typeof(c))/betamax
        else
            W[j,j] = 1.0/c
        end
    end
end

function init_optim_components(estimate::DiffEstimate, weights_type::AbstractNormWeights)
    @unpack uinit = estimate
    @unpack Nx = estimate.params.mesh
    Gc = zeros(eltype(uinit), Nx)
    A = zeros(eltype(uinit), 2*Nx, Nx)
    b = zeros(eltype(uinit), 2*Nx)
    W = zeros(eltype(uinit), Nx, Nx)

    Gflux!(CenteredG(), Gc, estimate)
    fill_A!(A, estimate)
    fill_b!(b, estimate)
    fill_W!(weights_type, W, estimate)

    Gc, A, b, W
end