struct Rusanov <: SpaceScheme end

struct RusanovCache <: Cache
    sL::Int
    sR::Int
    A
    
    function RusanovCache()
        new(1, 1, 0.0)
    end
end

function num_flux(::Rusanov, )
    uL, uR = u[1,:], u[2,:]
    A = CFL_cond(equation, u)
    (flux(equation, uL) .+ flux(equation, uR)) / 2 .- A / 2 * (uR .- uL)
end