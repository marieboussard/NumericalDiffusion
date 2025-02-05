struct MixedScheme{T<:AbstractFloat,U} <: SpaceScheme
    # (1-alpha)*f(scheme1) + alpha*f(scheme2)
    CFL_factor::T
    scheme1::SpaceScheme
    scheme2::SpaceScheme
    mixingWeights::AbstractArray{U} #where U
end
get_sL(ms::MixedScheme) = max(get_sL(ms.scheme1), get_sL(ms.scheme2))
get_sR(ms::MixedScheme) = max(get_sR(ms.scheme1), get_sR(ms.scheme2))
get_name(ms::MixedScheme) = "Mixed btw " * get_name(ms.scheme1) * " and " * get_name(ms.scheme2)


function vecNumFlux(::NullSource, ms::MixedScheme, equation::Equation, v; kwargs...)
    Nx, p = size(v, 1), get_unknowns_number(equation)
    sL, sR = get_sL(scheme), get_sR(scheme)
    alpha = ms.mixingWeights
    numericalFluxMat = zeros(eltype(v), Nx + 1, p)
    for i ∈ 2:Nx+1
        numericalFluxMat[i, :] .= (1 - alpha[i]) * numFlux(ms.scheme1, equation, extract_data_stencil(v, i - 1, sL, sR); kwargs...) + alpha[i] * numFlux(ms.scheme2, equation, extract_data_stencil(v, i - 1, sL, sR); kwargs...)
    end
    numericalFluxMat[1, :] = numericalFluxMat[end, :]
    numericalFluxMat
end


# function giveNumFlux(::NullSource, ms::MixedScheme{T,U}, equation::Equation, v; kwargs...) where {T} where {U}
#     #println("Using Mixed Scheme giveNumFlux function")
#     alpha = ms.mixingWeights
#     #T = eltype(alpha)
#     Nx, p = length(v[:, 1]), get_unknowns_number(equation)
#     numericalFluxMat = zeros(U, (Nx + 1, p))
#     for i ∈ 2:Nx
#         numericalFluxMat[i, :] = (1 - alpha[i]) * giveNumFlux(ms.scheme1, equation, v[i-1, :], v[i, :]) + alpha[i] * giveNumFlux(ms.scheme2, equation, v[i-1, :], v[i, :])
#     end
#     #@code_warntype giveNumFlux(ms.scheme1, equation, v[end,:], v[1,:])
#     numericalFluxMat[1, :] = (1 - alpha[1]) * giveNumFlux(ms.scheme1, equation, v[end, :], v[1, :]) + alpha[1] * giveNumFlux(ms.scheme2, equation, v[end, :], v[1, :])
#     numericalFluxMat[end, :] = numericalFluxMat[1, :]
#     numericalFluxMat
# end