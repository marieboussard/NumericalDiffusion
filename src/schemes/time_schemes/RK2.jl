struct RK2 <: TimeScheme
end

get_sL(::RK2, sp::SpaceScheme) = 2*get_sL(sp::SpaceScheme)
get_sR(::RK2, sp::SpaceScheme) = 2*get_sR(sp::SpaceScheme)

# get_sL(::RungeKutta, sp::SpaceScheme) = 1*get_sL(sp::SpaceScheme)
# get_sR(::RungeKutta, sp::SpaceScheme) = 1*get_sR(sp::SpaceScheme)

get_name(::RK2) = "RK2"

function giveNumFlux(::RK2, spaceScheme::SpaceScheme, equation::Equation, v...; dt, domain::Domain)
    vb = next_timestep(FVScheme(Euler(), spaceScheme), v, dt, domain, equation)
    (NumFlux(spaceScheme, equation, v[1+sL:1+2*sL+sR]...) + NumFlux(spaceScheme, equation, vb[1+sL:1+2*sL+sR]...))*0.5
end

# function next_timestep(::RK2, v, dt, domain::Domain, equation::Equation, spaceScheme::SpaceScheme)
#     vb = next_timestep(Euler(), v, dt, domain, equation, spaceScheme)
#     vbb = next_timestep(Euler(), vb, dt, domain, equation, spaceScheme)
#     (vbb .+ v)/2
# end

function exactG(::RK2, spaceScheme::SpaceScheme, equation::Equation, u, dt, domain::Domain)
    ub = next_timestep(Euler(), u, dt, domain, equation, spaceScheme)

    G = exactG(spaceScheme, equation, u)
    Gb = exactG(spaceScheme, equation, ub)

    (G .+ Gb)/2
end