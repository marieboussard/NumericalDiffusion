abstract type AbstractEntNumFlux end
struct CenteredG <: AbstractEntNumFlux end

function Gflux!(::CenteredG, res::AbstractVector, estimate::DiffEstimate)
    @unpack uinit, entfun = estimate
    @unpack Nx = estimate.params.mesh
    Gcont = zeros(eltype(uinit), Nx)
    G!(entfun, uinit, Gcont)
    for j in 1:Nx
        res[j] = 0.5*(Gcont[j]+Gcont[mod1(j+1,Nx)])
    end
end