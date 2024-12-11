abstract type Limiter end
struct Minmod<:Limiter end

struct MUSCL{T <: AbstractFloat} <: FVMethod
    CFL_factor::T
    dx::T
    subMethod::FVMethod
    limiter::Limiter
end

MUSCL(CFL_factor::Float64, subMethod::FVMethod, limiter::Limiter, domain::Domain) = MUSCL(CFL_factor, domain.dx, subMethod, limiter)

MUSCL(CFL_factor::Float64, subMethod::DataType, limiter::Limiter, domain::Domain) = MUSCL(CFL_factor, domain.dx, subMethod(CFL_factor), limiter)

get_sL(::MUSCL) = 2
get_sR(::MUSCL) = 2

get_name(::MUSCL) = "MUSCL"

# First we define different limiters (and associated slope)
compute_slope(::Minmod, uL, uC, uR) = max(0, min(uC-uL, uR-uC)) + min(0, max(uC-uL, uR-uC))

# # Then we define the affine interpolation
# affine_interp(a, b, x0, x) = a*(x-x0) + b

function numFlux(mu::MUSCL, equation::Equation, uL, uC, uR, uRR)

    # Piecewise linear reconstruction
    σC = compute_slope(mu.limiter, uL, uC, uR)
    σR = compute_slope(mu.limiter, uC, uR, uRR)
    u_minus = uC + mu.dx/2 * σC
    u_plus = uR + mu.dx/2 * σR

    numFlux(mu.subMethod, equation, u_minus, u_plus)

end