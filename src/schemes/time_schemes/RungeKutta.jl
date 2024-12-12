struct RungeKutta <: TimeScheme
end

get_sL(::RungeKutta, sp::SpaceScheme) = 2*get_sL(sp::SpaceScheme)
get_sR(::RungeKutta, sp::SpaceScheme) = 2*get_sR(sp::SpaceScheme)

get_name(::RungeKutta) = "Runge Kutta"

function next_timestep(::RungeKutta, v, dt, domain::Domain, equation::Equation, spaceScheme::SpaceScheme)
    vb = next_timestep(Euler(), v, dt, domain, equation, spaceScheme)
    vbb = next_timestep(Euler(), vb, dt, domain, equation, spaceScheme)
    (vbb .+ v)/2
end