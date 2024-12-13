struct RungeKutta <: TimeScheme
end

get_sL(::RungeKutta, sp::SpaceScheme) = 2*get_sL(sp::SpaceScheme)
get_sR(::RungeKutta, sp::SpaceScheme) = 2*get_sR(sp::SpaceScheme)

# get_sL(::RungeKutta, sp::SpaceScheme) = 1*get_sL(sp::SpaceScheme)
# get_sR(::RungeKutta, sp::SpaceScheme) = 1*get_sR(sp::SpaceScheme)

get_name(::RungeKutta) = "Runge Kutta"

function next_timestep(::RungeKutta, v, dt, domain::Domain, equation::Equation, spaceScheme::SpaceScheme)
    vb = next_timestep(Euler(), v, dt, domain, equation, spaceScheme)
    vbb = next_timestep(Euler(), vb, dt, domain, equation, spaceScheme)
    (vbb .+ v)/2
end

function exactG(::RungeKutta, spaceScheme::SpaceScheme, equation::Equation, u, dt, domain::Domain)
    ub = next_timestep(Euler(), u, dt, domain, equation, spaceScheme)

    G = exactG(spaceScheme, equation, u)
    Gb = exactG(spaceScheme, equation, ub)

    (G .+ Gb)/2
end