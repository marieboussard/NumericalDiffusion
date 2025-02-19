# include("abstract_methods.jl")
# include("../equations/source_terms.jl")

struct HR{T <: AbstractFloat} <: SpaceScheme
    CFL_factor::T
    subScheme::SpaceScheme
end

HR(CFL_factor::Float64, subScheme::DataType) = HR(CFL_factor, subScheme(CFL_factor))

get_sL(::HR) = 1
get_sR(::HR) = 1

get_name(::HR) = "Hydrostatic Reconstruction"

# function hMinus(vL, xL, xR, zbSource::ZbSource)
#     zL, zR = zb(zbSource, xL), zb(zbSource, xR)
#     hminus = max(0, vL[1] + zL - max(zL, zR))
# end

hMinus(vL, zL, zR) = max(0, vL[1] .+ zL .- max(zL, zR))
hPlus(vR, zL, zR) = max(0, vR[1] .+ zR .- max(zL, zR))

# function hPlus(vR, xL, xR, zbSource::ZbSource)
#     zL, zR = zb(zbSource, xL), zb(zbSource, xR)
#     hplus = max(0, vR[1] + zR - max(zL, zR))
# end

function numFlux(hydro::HR, equation::Equation, v; z=[0.0, 0.0], kwargs...)
    vL, vR = v[1,:], v[2,:]
    zL, zR = z
    if vL[1] > 1e-10
        #hminus = hMinus(vL, xL, xR, equation.source)
        hminus = hMinus(vL, zL, zR)
        vminus = [hminus, hminus * vL[2] / vL[1]]
    else
        vminus = [0.0, 0.0]
    end
    if vR[1] > 1e-10
        #hplus = hPlus(vR, xL, xR, equation.source)
        hplus = hPlus(vR, zL, zR)
        vplus = [hplus, hplus * vR[2] / vR[1]]
    else
        vplus = [0.0, 0.0]
    end
    # @show vminus, vplus
    # @show vcat(vminus', vplus')
    numFlux(hydro.subScheme, equation, vcat(vminus', vplus'))
end

# giveNumFlux(hydro::HR, equation::Equation, vL, vR; kwargs...) = numFlux(hydro, equation::Equation, vL, vR; kwargs...)
# function sourceTerm(::Hydrostatic, zbSource::ZbSource, domain::Domain, v)
#     z, dx = domain.sourceVec, domain.dx
#     Nx = size(v)[1]
#     result = zero(v)
#     for i in 1:Nx
#         result[i,1] = 0.0
#         result[i,2] = 1 / dx * g / 2 * ((hMinus(v[mod1(i, Nx),:], x[mod1(i, Nx)], x[mod1(i + 1, Nx)], zbSource))^2 - (hPlus(v[mod1(i, Nx),:], x[mod1(i - 1, Nx)], x[mod1(i, Nx)], zbSource))^2)
#     end
#     result
#     #1 / dx .* [[0.0, g / 2 * ((hMinus(v[mod1(i, Nx)], x[mod1(i, Nx)], x[mod1(i + 1, Nx)], zbSource))^2 - (hPlus(v[mod1(i, Nx)], x[mod1(i - 1, Nx)], x[mod1(i, Nx)], zbSource))^2)] for i in 1:Nx]
# end

# function sourceTerm(::Hydrostatic, zbSource::ZbSource, domain::Domain, v)
#     z, dx = domain.sourceVec, domain.dx
#     Nx = size(v)[1]
#     result = zero(v)
#     for i in 1:Nx
#         result[i,1] = 0.0
#         result[i,2] = 1/dx*g/2* ((hMinus(v[mod1(i, Nx),:], z[mod1(i, Nx)], z[mod1(i + 1, Nx)]))^2 - (hPlus(v[mod1(i, Nx),:], z[mod1(i - 1, Nx)], z[mod1(i, Nx)]))^2)
#     end
#     result
# end

function sourceTerm(::SaintVenant, ::HR, domain::Domain, v; z=domain.sourceVec)
    dx = domain.dx
    Nx = size(v)[1]
    result = zero(v)
    for i in 1:Nx
        result[i,1] = 0.0
        result[i,2] = 1/dx*g/2* ((hMinus(v[mod1(i, Nx),:], z[mod1(i, Nx)], z[mod1(i + 1, Nx)]))^2 - (hPlus(v[mod1(i, Nx),:], z[mod1(i - 1, Nx)], z[mod1(i, Nx)]))^2)
    end
    result
end