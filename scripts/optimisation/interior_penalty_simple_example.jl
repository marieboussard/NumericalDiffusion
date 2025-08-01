struct InteriorPenaltyCache{wtype<:AbstractMatrix, gtype<:AbstractVector, atype<:AbstractMatrix, btype<:AbstractVector, ttype<:AbstractArray} <: Cache
    # Parameters
    Nx::Int
    eps::Real

    # Problem components
    W::wtype
    Gc::gtype
    A::atype
    b::btype

    # Temporary storage
    Wgamma::gtype
    WGc::gtype 
    T::ttype 
    Agamma::btype

    function InteriorPenaltyCache(Nx::Int, eps::Real, W::AbstractMatrix, Gc::AbstractVector, A::AbstractMatrix, b::AbstractVector)
        Wgamma = zero(Gc)
        WGc = zero(Gc)
        T = zeros(1,Nx)
        Agamma = zero(b)
        new{typeof(W), typeof(Gc), typeof(A), typeof(b), typeof(T)}(Nx, eps, W, Gc, A, b, Wgamma, WGc, T, Agamma)
    end

end

function J_interior_penalty(gamma::AbstractVector, cache::InteriorPenaltyCache)
    @unpack Nx, eps, W, Gc, A, b, Wgamma, WGc, T, Agamma = cache
    mul!(Wgamma, W, gamma)
    mul!(WGc, W, Gc)
    transpose!(T, Wgamma.-WGc)
    res = (T * (gamma .- Gc))[1]
    mul!(Agamma, A, gamma)
    for i in 1:size(A)[1]
        res += eps/(Agamma[i]-b[i])
    end
    res
end