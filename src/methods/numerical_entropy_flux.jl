abstract type NumEntFlux <: Method end

struct CenteredG <: NumEntFlux end

numFlux(::CenteredG, equation::Equation, uL, uR) = (get_G(equation, uL)+get_G(equation, uR))/2


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