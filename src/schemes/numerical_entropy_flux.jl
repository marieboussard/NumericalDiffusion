abstract type NumEntFlux <: Scheme end

struct CenteredG <: NumEntFlux end

get_sL(::CenteredG) = 1
get_sR(::CenteredG) = 1

numFlux(::CenteredG, equation::Equation, u; kwargs...) = (G(equation, u[begin]) + G(equation, u[end])) / 2


# function giveAnalyticG(::CenteredG, u)
#     Nx, p = length(v[:,1]), get_unknowns_number(equation)
#     numericalFluxMat = zeros(eltype(v), Nx+1, p)
#     for i ∈ 2:Nx
#         numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:])
#     end
#     numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:])
#     numericalFluxMat[end,:] = numericalFluxMat[1,:]
#     #@show typeof(numericalFluxMat)
#     numericalFluxMat
# end

# function giveNumFlux(::NullSource, method::FVMethod, equation::Equation, v; kwargs...)
#     #println("Calling regular giveNumFlux")
#     Nx, p = length(v[:,1]), get_unknowns_number(equation)
#     numericalFluxMat = zeros(eltype(v), Nx+1, p)
#     for i ∈ 2:Nx
#         numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:])
#     end
#     numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:])
#     numericalFluxMat[end,:] = numericalFluxMat[1,:]
#     #@show typeof(numericalFluxMat)
#     numericalFluxMat
# end