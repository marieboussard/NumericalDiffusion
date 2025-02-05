function checkInequalities(domain::Domain, equation::Equation, scheme::FVScheme, u_approx, gamma, dt_vec)
    Nx = domain.Nx
    z = manageSource(domain)
    u, up = u_approx[end-1], u_approx[end]
    dx, dt = domain.dx, dt_vec[end]
    sL, sR = get_sL(scheme), get_sR(scheme)
    numericalFluxMat = vecNumFlux(equation.source, scheme, equation, u; dt=dt, domain=domain)

    e1, e2 = zeros(eltype(gamma), Nx), zeros(eltype(gamma), Nx)

    for j in 1:Nx

        uh_j = u[j,:] - dt/dx * (numericalFluxMat[j+1,:] - flux(equation, u[j,:]))
        uh_jp = u[mod1(j+1,Nx),:] - dt/dx * (flux(equation, u[mod1(j+1,Nx),:]) - numericalFluxMat[j+1,:])

        # TO DO: also compute the modified source term

        # @show get_eta(equation, uh_j, z[j])[1]
        # @show get_eta(equation, u[j,:], z[j])[1]

        #e1[j] = get_eta(equation, uh_j, z[j])[1] - get_eta(equation, u[j,:], z[j])[1] + dt / dx * (gamma[j+1] - get_G(equation, u[j,:], z[j])[1])
        e1[j] = get_eta(equation, uh_j, 0.0)[1] - get_eta(equation, u[j,:], 0.0)[1] + dt / dx * (gamma[j+1] - get_G(equation, u[j,:], 0.0)[1])

        e2[j] = get_eta(equation, uh_jp, z[mod1(j+1,Nx)])[1] - get_eta(equation, u[mod1(j+1,Nx),:], z[mod1(j+1,Nx)])[1] + dt / dx * (get_G(equation, u[mod1(j+1,Nx),:], z[mod1(j+1,Nx)])[1] - gamma[j+1])

        # ut = compute_u_tilde(modifiedDataType, u, j, sL, sR)
        # zt = isnothing(domain.sourceVec) ? zero(ut) : compute_u_tilde(modifiedDataType, z, j, sL, sR)
        # uh = compute_u_hat(equation.source, ut, dx, dt, j, domain, equation, method; zt=zt)

        # e1[j] = get_eta(equation, uh[j,:], zt[j])[1] - get_eta(equation, ut[j,:], zt[j])[1] + dt / dx * (gamma[j+1] - get_G(equation, ut[j,:], zt[j])[1])

        # e2[j] = get_eta(equation, uh[mod1(j+1,Nx),:], zt[mod1(j+1,Nx)])[1] - get_eta(equation, ut[mod1(j+1,Nx),:], zt[mod1(j+1,Nx)])[1] + dt / dx * (get_G(equation, ut[j,:], zt[j])[1] - gamma[j+1])


        # e1 = [get_eta(equation, up[i,:], z[i])[1] - get_eta(equation, u[i,:], z[i])[1] + dt / dx * (gamma[i+1] - get_G(equation, u[i,:], z[i])[1]) for i in 1:length(u[:,1])]

        # e2 = [get_eta(equation, up[mod1(i+1,Nx),:], z[mod1(i+1,Nx)])[1] - get_eta(equation, u[mod1(i+1,Nx),:], z[mod1(i+1,Nx)])[1] + dt / dx * (get_G(equation, u[i,:], z[i])[1] - gamma[i+1]) for i in 1:Nx]
    end
    e1, e2
end

checkInequalities(solEnt::OptForEntropySol) = checkInequalities(solEnt.domain, solEnt.equation, solEnt.scheme, solEnt.u_approx, solEnt.Gopt, solEnt.dt_vec)