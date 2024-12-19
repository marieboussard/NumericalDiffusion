function extract_data_stencil(u, j, sL, sR)
    Nx = size(u)[1]
    if length(size(u)) > 1
        p = size(u)[2]
    else
        p = 1
    end
    u_short = zeros(eltype(u), (sL+sR,p))
    i = 1
    for k in j-sL+1:j+sR 
        u_short[i,:] = u[mod1(k, Nx),:]
        i+=1
    end
    u_short
end

function vecNumFlux(::NullSource, scheme::Scheme, equation::Equation, v; kwargs...)
    Nx, p = size(v, 1), get_unknowns_number(equation)
    sL, sR = get_sL(scheme), get_sR(scheme)
    numericalFluxMat = zeros(eltype(v), Nx+1, p)
    for i ∈ 2:Nx+1
        numericalFluxMat[i,:] .= numFlux(scheme, equation, extract_data_stencil(v, i-1, sL, sR); kwargs...)
    end
    numericalFluxMat[1,:] = numericalFluxMat[end,:]
    numericalFluxMat
end

numFlux(scheme::FVScheme, args...; kwargs...) = numFlux(scheme.timeScheme, scheme.spaceScheme, args...; kwargs...)

function vecNumFlux(::ZbSource, scheme::Scheme, equation::Equation, v; domain::Domain, kwargs...)
    Nx, p = size(v, 1), get_unknowns_number(equation)
    sL, sR = get_sL(scheme), get_sR(scheme)
    numericalFluxMat = zeros(eltype(v), Nx+1, p)
    for i ∈ 2:Nx+1
        numericalFluxMat[i,:] = numFlux(scheme, equation, extract_data_stencil(v, i-1, sL, sR); z = extract_data_stencil(domain.sourceVec, i-1, sL, sR), kwargs...)
    end
    numericalFluxMat[1,:] = numericalFluxMat[end,:]
    numericalFluxMat
end