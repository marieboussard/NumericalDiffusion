function extract_data_stencil(equation::Equation, u, j, sL, sR)
    Nx, p = length(u[:,1]), get_unknowns_number(equation)
    u_short = zeros(eltype(u), (sL+sR,p))
    i = 1
    for k in j-sL+1:j+sR 
        u_short[i,:] = u[mod1(k, Nx),:]
        i+=1
    end
    u_short
end

function giveNumFlux(::NullSource, method::Method, equation::Equation, v; kwargs...)
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    sL, sR = get_sL(method), get_sR(method)
    numericalFluxMat = zeros(eltype(v), Nx+1, p)
    for i ∈ 2:Nx+1
        #@show (extract_data_stencil(equation, v, i, sL, sR)...)
        # if i!=1
        #     @show giveNumFlux(method, equation, v[i-1,:], v[i,:])
        # end
        #@show giveNumFlux(method, equation, extract_data_stencil(equation, v, i-1, sL, sR)...)
        numericalFluxMat[i,:] .= giveNumFlux(method, equation, extract_data_stencil(equation, v, i-1, sL, sR)...)
    end
    # numericalFluxMat[1,:] = giveNumFlux(method, equation, v[end,:], v[1,:])
    # numericalFluxMat[end,:] = numericalFluxMat[1,:]
    numericalFluxMat[1,:] = numericalFluxMat[end,:]

    numericalFluxMat
end







# function giveNumFlux(::NullSource, method::Method, equation::Equation, v; kwargs...)
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