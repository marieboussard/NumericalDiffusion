struct MixedMethod <: FVMethod
    # (1-alpha)*f(method1) + alpha*f(method2)
    CFL_factor::Float64
    method1::FVMethod
    method2::FVMethod
    mixingWeights
end
get_sL(mm::MixedMethod) = max(get_sL(mm.method1), get_sL(mm.method2))
get_sR(mm::MixedMethod) = max(get_sR(mm.method1), get_sR(mm.method2))
get_name(mm::MixedMethod) = "Mixed btw "*get_name(mm.method1)*" and "*get_name(mm.method2)

function giveNumFlux(::NullSource, mm::MixedMethod, equation::Equation, u; kwargs...)
    alpha = mm.mixingWeights
    Nx, p = length(v[:,1]), get_unknowns_number(equation)
    numericalFluxMat = zeros(Nx+1, p)
    for i ∈ 2:Nx
        numericalFluxMat[i,:] = (1-alpha[i])*giveNumFlux(mm.method1, equation, v[i-1,:], v[i,:]) 
        + alpha[i]*giveNumFlux(mm.method2, equation, v[i-1,:], v[i,:])
    end
    numericalFluxMat[1,:] = (1-alpha[1])*giveNumFlux(mm.method1, equation, v[end,:], v[1,:])
    + alpha[1]*giveNumFlux(mm.method2, equation, v[end,:], v[1,:])
    numericalFluxMat[end,:] = numericalFluxMat[1,:]
    numericalFluxMat
end