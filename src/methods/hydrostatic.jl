include("abstract_methods.jl")
include("../equations/source_terms.jl")

struct Hydrostatic <: FVMethod
    CFL_factor::Float64
    subMethod::FVMethod
end

get_sL(::Hydrostatic) = 1
get_sR(::Hydrostatic) = 1

function hMinus(vL, xL, xR, zbSource::ZbSource)
    zL, zR = zb(zbSource, xL), zb(zbSource, xR)
    hminus = max(0, vL[1] + zL - max(zL, zR))
end

function hPlus(vR, xL, xR, zbSource::ZbSource)
    zL, zR = zb(zbSource, xL), zb(zbSource, xR)
    hplus = max(0, vR[1] + zR - max(zL, zR))
end

function numFlux(hydro::Hydrostatic, equation::Equation, vL, vR; xL=0, xR=0)
    hminus, hplus = hMinus(vL, xL, xR, equation.source), hPlus(vR, xL, xR, equation.source)
    vminus, vplus = [hminus, hminus * vL[2] / vL[1]], [hplus, hplus * vR[2] / vR[1]]
    numFlux(hydro.subMethod, equation, vminus, vplus)
end

giveNumFlux(hydro::Hydrostatic, equation::Equation, vL, vR; kwargs...) = numFlux(hydro, equation::Equation, vL, vR; kwargs...)
