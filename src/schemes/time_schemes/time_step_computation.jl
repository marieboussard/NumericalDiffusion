function next_timestep(scheme::FVScheme, v, dt, domain::Domain, equation::Equation)
    numericalFluxMat = giveNumFlux(NullSource(), scheme, equation, collect(v); dt=dt, domain=domain)
    v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
end