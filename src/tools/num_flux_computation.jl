function giveNumFlux(::NullSource, method::Method, equation::Equation, v; kwargs...)
    #println("Calling regular giveNumFlux")
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(eltype(v), Nx+1, p)
    for i ∈ 2:Nx
        numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    #@show typeof(numericalFluxMat)
    numericalFluxMat
end

function giveNumFlux(::ZbSource, method::Method, equation::Equation, v; domain::Domain)
    #println("Entering giveNumFlux")
    #@show eltype(v)
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(eltype(v), Nx+1, p)
    for i ∈ 2:Nx
        numericalFluxMat[i,:] = giveNumFlux(method, equation, v[i-1,:], v[i,:]; zL=domain.sourceVec[i-1], zR=domain.sourceVec[i])
    end
    numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:]; zL=domain.sourceVec[end], zR=domain.sourceVec[1])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    #@show typeof(numericalFluxMat)
    numericalFluxMat
end