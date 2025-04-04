abstract type Limiter end
struct Minmod <: Limiter end
struct Superbee <: Limiter end

struct MUSCL{T<:AbstractFloat} <: SpaceScheme
    CFL_factor::T
    subScheme::SpaceScheme
    limiter::Limiter
end

# MUSCL(CFL_factor::Float64, subScheme::SpaceScheme, limiter::Limiter) = MUSCL(CFL_factor, subScheme, limiter)

MUSCL(CFL_factor::Float64, subScheme::DataType, limiter::Limiter, domain::Domain) = MUSCL(CFL_factor, domain.dx, subScheme(CFL_factor), limiter)

get_sL(::MUSCL) = 2
get_sR(::MUSCL) = 2

get_name(::MUSCL) = "MUSCL"

# First we define different limiters (and associated slope)
minmod(x, y) = max(0, min(x, y)) + min(0, max(x, y))
maxmod(x, y) = max(0, max(x, y)) + min(0, min(x, y))
s1(uL, uC, uR) = minmod(uR - uC, 2 * (uC - uL))
s2(uL, uC, uR) = minmod(2 * (uR - uC), uC - uL)
#compute_slope(::Minmod, uL, uC, uR) = max(0, min(uC-uL, uR-uC)) + min(0, max(uC-uL, uR-uC))
compute_slope(::Minmod, uL, uC, uR) = minmod(uC - uL, uR - uC)
compute_slope(::Superbee, uL, uC, uR) = maxmod(s1(uL, uC, uR), s2(uL, uC, uR))

# # Then we define the affine interpolation
# affine_interp(a, b, x0, x) = a*(x-x0) + b

function numFlux(mu::MUSCL, equation::Equation, u; kwargs...)
    #println("computing muscl flux")

    #@show uL, uC, uR, uRR

    uL, uC, uR, uRR = u

    # Piecewise linear reconstruction
    σC = compute_slope(mu.limiter, uL, uC, uR)
    σR = compute_slope(mu.limiter, uC, uR, uRR)
    #@show mu.dx
    uC_plus = uC + 1 / 2 * σC
    uR_minus = uR - 1 / 2 * σR

    #@show u_plus, u_minus

    #@show mu.subScheme

    # @show uC, uR
    # @show uC_plus, uR_minus

    #@show numFlux(mu.subScheme, equation, uC, uR)
    numFlux(mu.subScheme, equation, [uC_plus, uR_minus]; kwargs...)

end