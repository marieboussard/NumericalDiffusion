struct MixedMethod <: FVMethod
    CFL_factor::Float64
    method1::FVMethod
    method2::FVMethod
    mixingWeights
end
get_sL(mm::MixedMethod) = max(get_sL(mm.method1), get_sL(mm.method2))
get_sR(mm::MixedMethod) = max(get_sR(mm.method1), get_sR(mm.method2))
get_name(mm::MixedMethod) = "Mixed btw "*get_name(mm.method1)*" and "*get_name(mm.method2)

function numFlux(mm:MixedMethod, equation::Equation, uL, uR)
end