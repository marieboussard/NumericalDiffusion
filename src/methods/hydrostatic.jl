# include("abstract_methods.jl")
# include("../equations/source_terms.jl")

struct Hydrostatic <: FVMethod
    CFL_factor::Float64
    subMethod::FVMethod
end

createHydrostatic(CFL_factor::Float64, subMethod::DataType) = Hydrostatic(CFL_factor, subMethod(CFL_factor))

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
function sourceTerm(::Hydrostatic, zbSource::ZbSource, domain::Domain, v)
    x, dx, Nx = domain.x, domain.dx, domain.Nx
    result = zero(v)
    for i in 1:Nx
        result[i,1] = 0.0
        result[i,2] = 1 / dx * g / 2 * ((hMinus(v[mod1(i, Nx),:], x[mod1(i, Nx)], x[mod1(i + 1, Nx)], zbSource))^2 - (hPlus(v[mod1(i, Nx),:], x[mod1(i - 1, Nx)], x[mod1(i, Nx)], zbSource))^2)
    end
    result
    #1 / dx .* [[0.0, g / 2 * ((hMinus(v[mod1(i, Nx)], x[mod1(i, Nx)], x[mod1(i + 1, Nx)], zbSource))^2 - (hPlus(v[mod1(i, Nx)], x[mod1(i - 1, Nx)], x[mod1(i, Nx)], zbSource))^2)] for i in 1:Nx]
end