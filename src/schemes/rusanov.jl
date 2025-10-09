mutable struct RusanovCache{utype<:AbstractArray} <: SpaceCache
    A::Float64
    fcont::utype
    absDfcont::utype
    # function RusanovCache(u::AbstractVector, jstart::Int=1, jend::Int=length(u))
    #     # new(zero(Float64), zero(u[jstart:jend]))
    #     fcont = zeros(eltype(u), jend-jstart+1)
    #     absDfcont = zero(fcont)
    #     new(zero(Float64), fcont, absDfcont)
    # end
    function RusanovCache(u::AbstractArray)
        fcont = zero(u)
        absDfcont = zero(fcont)
        new{typeof(fcont)}(zero(Float64), fcont, absDfcont)
    end
end



struct Rusanov <: SpaceScheme end

get_sL(::Rusanov) = 1
get_sR(::Rusanov) = 1

get_name(::Rusanov) = "Rusanov"


function update_cache!(rcache::RusanovCache, u::AbstractArray, eqfun::AbstractEquationFun, jstart::Int=1, jend::Int=length(u))
    @unpack fcont, absDfcont = rcache
    flux!(eqfun, view(u, jstart:jend), fcont)
    Dflux!(eqfun, view(u, jstart:jend), absDfcont)
    abs!(absDfcont, absDfcont)
end

# const g = 9.8
# const treshold = 1e-6 # treshold for dealing with dry water states (h=0)

function update_cache!(rcache::RusanovCache, u::AbstractArray, ::SaintVenant, jstart::Int=1, jend::Int=size(u)[1])
    @unpack absDfcont = rcache 
    for j in jstart:jend
        if u[j,1] > treshold
            absDfcont[j] = abs(u[j,2]/ u[j,1]) + sqrt(g*u[j,1])
        else
            absDfcont[j] = 0.0
        end
    end
end


# Diffusion coefficient A(uL, uR) (depends on the equation)
function ARusanov!(::OneD, ::EquationType, ::AbstractEquationFun, rcache::RusanovCache, j::Int, Nx::Int)
    @unpack absDfcont = rcache
    rcache.A = absDfcont[j]
    rcache.A = max(rcache.A, absDfcont[mod1(j+1,Nx)])
end

# function ARusanov!(::OneD, ::System, ::SaintVenant, rcache::RusanovCache, j::Int, Nx::Int)
#     @unpack absDfcont = rcache
#     rcache.A = absDfcont[j]
#     rcache.A = max(rcache.A, absDfcont[mod1(j+1,Nx)])
# end


function numflux!(::Rusanov, rcache::RusanovCache, equation::Equation, u::AbstractArray, fnum::AbstractArray, ju::Int, Nx=size(u)[1], jf::Int=ju)
    @unpack fcont = rcache
    ARusanov!(equation.dim, equation.eqtype, equation.funcs, rcache, ju, Nx)
    for r in 1:equation.p
        fnum[jf,r] = (fcont[ju, r] + fcont[mod1(ju+1,Nx),r]) *0.5 - rcache.A.*0.5 * (u[mod1(ju+1,Nx),r] - u[ju,r])
    end
end





# function numflux!(::Rusanov, integrator::Integrator, j::Int, args...)
#     @unpack equation, params, cache, space_cache, fnum, fcont, uprev = integrator
#     @unpack Nx = params.mesh
#     CFL_local!(equation.dim, equation.eqtype, equation.funcs, integrator, j)
#     for r in 1:equation.p
#         # fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
#         fnum[j,r] = (fcont[j, r] + fcont[mod1(j+1,Nx),r]) *0.5 - space_cache.cfl_loc.*0.5 * (uprev[mod1(j+1,Nx),r] - uprev[j,r])
#     end
# end

# function numflux!(::Rusanov, j::Int, params::Parameters, equation::Equation, cache::Cache, space_cache::SpaceCache, fnum::AbstractArray, fcont::AbstractArray, u::AbstractArray, i::Int=j)
#     @unpack Nx = params.mesh
#     CFL_local!(equation.dim, equation.eqtype, equation.funcs, j, params, cache, space_cache)
#     for r in 1:equation.p
#         fnum[i,r] = (fcont[j, r] + fcont[mod1(j+1,Nx),r]) *0.5 - space_cache.cfl_loc.*0.5 * (u[mod1(j+1,Nx),r] - u[j,r])
#     end
# end

# # function numflux!(::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, subcache::scacheType, integrator::Integrator, j::Int, args...)
# #     @unpack equation, params, fnum = integrator
# #     @unpack Nx = params.mesh
# #     subcache.cfl_loc = max(abs(uL[2]/uL[1]) + sqrt(g*uL[1]), abs(uR[2]/uR[1]) + sqrt(g*uR[1]))
# #     for r in 1:equation.p
# #         fnum[j,r] = (fL[r] + fR[r]) *0.5 - subcache.cfl_loc.*0.5 * (uR[r] - uL[r])
# #     end
# # end

# function numflux!(::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, fnum::AbstractVector, subcache::SpaceCache, equation::Equation, args...)
#     # This version is suitable for Saint-Venant only
#     subcache.cfl_loc = max(abs(uL[2]/uL[1]) + sqrt(g*uL[1]), abs(uR[2]/uR[1]) + sqrt(g*uR[1]))
#     for r in 1:equation.p
#         fnum[r] = (fL[r] + fR[r]) *0.5 - subcache.cfl_loc.*0.5 * (uR[r] - uL[r])
#     end
# end

# numflux!(scheme::Rusanov, uL::AbstractVector, uR::AbstractVector,  fL::AbstractVector, fR::AbstractVector, subcache::SpaceCache, integrator::Integrator, j::Int, args...) = numflux!(scheme, uL, uR, fL, fR, view(integrator.fnum,j,:), subcache, integrator.equation, args...)