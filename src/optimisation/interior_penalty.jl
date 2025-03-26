# struct InteriorPenCache <: Cache

#     # Problem components
#     W::wtype
#     Gc::gtype
#     A::atype
#     b::btype

#     # Temporary calculus
#     Agamma::btype

#     # Parameters
#     eps::Float64
# end

struct InteriorPenCache{wtype<:AbstractMatrix, gtype<:AbstractVector} <: Cache
    
    # Parameters
    Nx::Int
    dt::Float64
    dx::Float64
    eps::Float64
    
    # Problem components
    W::wtype
    Gc::gtype
    l::gtype
    L::gtype
    
end

function JIntPen(gamma::AbstractVector, cache::InteriorPenCache)
    res = zero(eltype(gamma))
    @unpack Nx, dt, dx, W, Gc, eps, L, l = cache
    for j in 1:Nx
        # @show dt/dx*(gamma[j]-gamma[mod1(j-1,Nx)])-L[j]
        # @show dt/dx*(gamma[mod1(j-1,Nx)]-gamma[j])+l[j]
        res += W[j]*norm(gamma .- Gc) - eps*(1.0/(dt/dx*(gamma[j]-gamma[mod1(j-1,Nx)])-L[j]) + 1.0/(dt/dx*(gamma[mod1(j-1,Nx)]-gamma[j])+l[j]))
    end
    res
end