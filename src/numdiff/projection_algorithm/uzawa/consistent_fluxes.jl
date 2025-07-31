abstract type AbstractEntNumFlux end
struct CenteredG <: AbstractEntNumFlux end
struct ConvexCombinationG{T<:Real} <: AbstractEntNumFlux
theta::T
end

function Gflux!(::CenteredG, res::AbstractVector, Nx::Int, uinit::AbstractArray, entfun::AbstractEntropyFun)
    Gcont = zeros(eltype(uinit), Nx)
    G!(entfun, uinit, Gcont)
    for j in 1:Nx
        res[j] = 0.5*(Gcont[j]+Gcont[mod1(j+1,Nx)])
    end
end

Gflux!(ent_num_flux::AbstractEntNumFlux, res::AbstractVector, estimate::DiffEstimate) = Gflux!(ent_num_flux, res, estimate.params.mesh.Nx, estimate.uinit, estimate.entfun)

# function Gflux!(::CenteredG, res::AbstractVector, estimate::DiffEstimate)
#     @unpack uinit, entfun = estimate
#     @unpack Nx = estimate.params.mesh
#     Gcont = zeros(eltype(uinit), Nx)
#     G!(entfun, uinit, Gcont)
#     for j in 1:Nx
#         res[j] = 0.5*(Gcont[j]+Gcont[mod1(j+1,Nx)])
#     end
# end

function Gflux!(ccg::ConvexCombinationG, res::AbstractVector, Nx::Int, uinit::AbstractArray, entfun::AbstractEntropyFun)
    Gcont = zeros(eltype(uinit), Nx)
    G!(entfun, uinit, Gcont)
    for j in 1:Nx
        res[j] = ccg.theta*Gcont[j] + (1-ccg.theta) * Gcont[mod1(j+1,Nx)]
    end
end