abstract type Limiter end
struct Minmod<:Limiter end

struct MUSCL{T <: AbstractFloat} <: SpaceScheme
    CFL_factor::T
    dx::T
    subScheme::SpaceScheme
    limiter::Limiter
end

MUSCL(CFL_factor::Float64, subScheme::SpaceScheme, limiter::Limiter, domain::Domain) = MUSCL(CFL_factor, domain.dx, subScheme, limiter)

MUSCL(CFL_factor::Float64, subScheme::DataType, limiter::Limiter, domain::Domain) = MUSCL(CFL_factor, domain.dx, subScheme(CFL_factor), limiter)

get_sL(::MUSCL) = 2
get_sR(::MUSCL) = 2

get_name(::MUSCL) = "MUSCL"

# First we define different limiters (and associated slope)
compute_slope(::Minmod, uL, uC, uR) = max(0, min(uC-uL, uR-uC)) + min(0, max(uC-uL, uR-uC))

# # Then we define the affine interpolation
# affine_interp(a, b, x0, x) = a*(x-x0) + b

function numFlux(mu::MUSCL, equation::Equation, uL, uC, uR, uRR)

    #@show uL, uC, uR, uRR

    # Piecewise linear reconstruction
    σC = compute_slope(mu.limiter, uL, uC, uR)
    σR = compute_slope(mu.limiter, uC, uR, uRR)
    u_plus = uC + mu.dx/2 * σC
    u_minus = uR - mu.dx/2 * σR

    #@show u_plus, u_minus

    numFlux(mu.subScheme, equation, u_plus, u_minus)

end