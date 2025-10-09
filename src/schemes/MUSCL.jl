# Limiters
abstract type Limiter end
struct Minmod <: Limiter end
struct Superbee <: Limiter end

# Different limiters (and associated slope)
minmod(x::Real, y::Real) = max(0, min(x, y)) + min(0, max(x, y))
maxmod(x::Real, y::Real) = max(0, max(x, y)) + min(0, min(x, y))
s1(uL::Real, uC::Real, uR::Real) = minmod(uR - uC, 2 * (uC - uL))
s2(uL::Real, uC::Real, uR::Real) = minmod(2 * (uR - uC), uC - uL)
#compute_slope(::Minmod, uL, uC, uR) = max(0, min(uC-uL, uR-uC)) + min(0, max(uC-uL, uR-uC))
compute_slope(::Minmod, uL::Real, uC::Real, uR::Real) = minmod(uC - uL, uR - uC)
compute_slope(::Superbee, uL::Real, uC::Real, uR::Real) = maxmod(s1(uL, uC, uR), s2(uL, uC, uR))

get_name(::Minmod) = "minmod"
get_name(::Superbee) = "superbee"

struct MUSCL{stype<:SpaceScheme, ltype<:Limiter} <: SpaceScheme
    subscheme::stype
    limiter::ltype
end

get_sL(::MUSCL) = 2
get_sR(::MUSCL) = 2

get_name(ms::MUSCL) = "MUSCL"*" ("*get_name(ms.limiter)*")"


# CACHE

mutable struct MUSCLCache{sctype<:SpaceCache, ltype<:Limiter, slptype<:AbstractVector, utype<:AbstractVector} <: SpaceCache 

    subcache::sctype
    limiter::ltype
    slope::slptype
    uLR::utype
    dx::Float64

    function MUSCLCache(ms::MUSCL, u::AbstractVector, dim::EquationDim, params::Parameters, args...)
        slope = zero(u)
        uLR = zeros(eltype(u), 2)
        subcache = init_cache(ms.subscheme, uLR, dim, params, args...)
        new{typeof(subcache), typeof(ms.limiter), typeof(slope), typeof(uLR)}(subcache, ms.limiter, slope, uLR, params.mesh.dx)
    end

end

function update_cache!(mcache::MUSCLCache, u::AbstractArray, ::AbstractEquationFun, jstart::Int=1, jend::Int=length(u))

    @unpack limiter, slope = mcache
    Nx = length(u)

    # Compute the slopes with the limiter

    for j in jstart:jend
        slope[j] = compute_slope(limiter, u[mod1(j-1,Nx)], u[j], u[mod1(j+1, Nx)])
    end

end


# NUMERICAL FLUX

function numflux!(scheme::MUSCL, mcache::MUSCLCache, equation::Equation, u::AbstractVector, fnum::AbstractVector, ju::Int, Nx=length(u), jf::Int=ju)

    @unpack subscheme = scheme
    @unpack slope, uLR, subcache, dx = mcache

    # Compute the boundary extrapolated values

    uLR[1] = u[ju] + slope[ju] / 2
    uLR[2] = u[mod1(ju + 1, Nx)] - slope[mod1(ju + 1, Nx)] / 2
    

    # Update low order cache

    update_cache!(subcache, uLR, equation.funcs)

    # Compute low order flux

    numflux!(subscheme, subcache, equation, uLR, fnum, 1, Nx, jf)

end