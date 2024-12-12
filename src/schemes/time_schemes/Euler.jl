struct Euler<:TimeScheme
end

get_sL(::Euler, sp::SpaceScheme) = get_sL(sp::SpaceScheme)
get_sR(::Euler, sp::SpaceScheme) = get_sR(sp::SpaceScheme)

function next_timestep(::Euler, v, dt, domain::Domain, equation::Equation, spaceScheme::SpaceScheme)
    numericalFluxMat = giveNumFlux(NullSource(), spaceScheme, equation, v)
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
end