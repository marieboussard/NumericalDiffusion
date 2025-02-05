struct Euler<:TimeScheme
end

get_name(::Euler) = "Euler"

get_sL(::Euler, sp::SpaceScheme) = get_sL(sp::SpaceScheme)
get_sR(::Euler, sp::SpaceScheme) = get_sR(sp::SpaceScheme)

numFlux(::Euler, spaceScheme::SpaceScheme, args...; kwargs...) = numFlux(spaceScheme::SpaceScheme, args...; kwargs...)

# function next_timestep(::Euler, v, dt, domain::Domain, equation::Equation, spaceScheme::SpaceScheme)
#     #println("Calling Euler for scheme"*get_name(spaceScheme))
#     numericalFluxMat = giveNumFlux(NullSource(), spaceScheme, equation, v)
#     #display(plot!(domain.interfaces, numericalFluxMat, label=get_name(spaceScheme)))

#     #vp = v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])

#     #display(plot!(domain.x, vp, label=get_name(spaceScheme)))

#    # @show v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
#     v - dt / domain.dx * (numericalFluxMat[2:end,:] - numericalFluxMat[1:end-1,:])
# end

exactG(::Euler, spaceScheme::SpaceScheme, equation::Equation, u, args...) = exactG(spaceScheme, equation, u)