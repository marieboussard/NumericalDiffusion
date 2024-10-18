include("../src/include_file.jl")

# 1 # Solving Burgers equation
xmin, xmax, Nx, t0, Tf = -2, 2, 10, 0, 0.4
CFL_factor = 0.5
domain = createInterval(xmin, xmax, Nx, t0, Tf)
equation = burgers()
method = Roe(CFL_factor)
#method = Rusanov(CFL_factor)
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
T_reached = 0.0
searchWhenAlreadyEntropic = true

u0_temp = u0
for k in 1:3
    global u0_temp
    global T_reached
    global solEnt
    global muopt
    global solEntRus

    dx = domain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u0_temp)
    domain.Tf = dt # Tf is modified so that only one timestep is performed
    T_reached += dt
    solEnt = optimize_for_entropy(u0_temp, domain, equation, method)
    solEntRus = optimize_for_entropy(u0_temp, domain, equation, Rusanov(CFL_factor))

    # AG_init = vcat(zero(solEnt.Gopt), solEnt.Gopt)
    # @show AG_init = reshape(AG_init, size(AG_init)[1])
    # _p = 0.0
    # cons(res, AG, p) = (res .= [constraint_AG(AG, solEnt)])
    # @show constraint_AG(AG_init, solEnt)

    if maximum(solEnt.Dopt) > 1e-8 || searchWhenAlreadyEntropic==true
        println("Modifiying the scheme to compensate positive numerical diffusion")

        solmuG = find_optimal_muG(solEnt);
        muopt, Gopt_mod = solmuG.mu, solmuG.G

        muG = vcat(muopt, Gopt_mod)
        muG = reshape(muG, size(muG)[1])
        @show constraint_muG(muG[1:end-1], solEnt)

        # # First constraint: Negative diffusion
        # f_delta = f_delta_from_A(Aopt, solEnt)
        up_mod = modified_scheme_step_mu(muopt, solEnt)
        # lambd = dt/dx
        # up = solEnt.u_approx[end]
        # bu_delta = zeros(Nx)
        # bg_delta = zeros(Nx)
        # cons=0.0
        # for j in 1:Nx
        #     bu_delta[j] = get_eta(equation, up_mod[j]) - get_eta(equation, up[j])
        #     bg_delta[j] = lambd*((Gopt_mod[j+1] - solEnt.Gopt[j+1]) - (Gopt_mod[j] - solEnt.Gopt[j]))
        #     cons += min(0, min(0,-solEnt.Dopt[j])-bu_delta[j]-bg_delta[j])^2
        # end
        # @show cons

        # z = isnothing(domain.sourceVec) ? zeros(Nx,1) : domain.sourceVec
        # plot(domain.x, [get_eta(equation, solEnt.u_approx[end-1][i])[1] for i in 1:length(solEnt.u_approx[end-1])])
        # display(plot!(domain.x, [get_eta(equation, solEnt.u_approx[end-1][i,:]; z=z[i])[1] for i in 1:length(solEnt.u_approx[end-1][:,1])]))

        # plot(domain.dx, dt / dx * (Gopt_mod[2:end] - Gopt_mod[1:end-1]))
        # display(plot!(domain.dx, ))

        #@show typeof([get_eta(equation, up_mod[j])[1].-get_eta(equation, solEnt.u_approx[end-1])[1]  for j in 1:length(up_mod)])
        #@show typeof(dt / dx * (Gopt_mod[2:end] - Gopt_mod[1:end-1]))

        # manual_D = [get_eta(equation, up_mod[j])[1]-get_eta(equation, solEnt.u_approx[end-1])[1]  for j in 1:length(up_mod)].+dt / dx * (Gopt_mod[2:end] .- Gopt_mod[1:end-1])


        # plot(domain.x, bu_delta, label="bu_delta")
        # plot!(domain.x, bg_delta, label="bg_delta")
        # display(plot!(domain.x, solEnt.Dopt, label="Dopt"))

        # plot(domain.x, bu_delta.+bg_delta, label="B")
        # plot!(domain.x, bu_delta.+bg_delta.+solEnt.Dopt, label="B+Dopt")
        # display(plot!(domain.x, -solEnt.Dopt, label="-Dopt"))

        u0_temp = up_mod

        #u0_temp = modified_scheme_step(Aopt, solEnt)
        # plot(domain.x, u0_temp, label="u0_temp")
        # display(plot!(domain.x, up_mod, label="up mod"))

        # m_delta, M_delta = compute_modified_bounds_f_delta(f_delta_from_A(Aopt, solEnt), solEnt)

        # plot(domain.interfaces, m_delta, label="m delta")
        # plot!(domain.interfaces, Gopt_mod, label="Gopt mod")
        # plot!(domain.interfaces, M_delta, label="M delta")
        # display(title!("Bounds at t="*string(round(T_reached, sigdigits=3))))

        plot(domain.interfaces, Gopt_mod[1:end-1], label="Gopt mod")
        plot!(domain.interfaces, solEnt.Gopt, label="Gopt"*get_name(method))
        plot!(domain.interfaces, solEntRus.Gopt, label="Gopt Rusanov")
        display(title!("Num entropy flux at t="*string(round(T_reached, sigdigits=3))))

        Dopt_mod = diffusion(solEnt.u_approx[end-1], up_mod, Gopt_mod[1:end-1], domain.dx, solEnt.dt_vec[end], equation, domain)
        
        # plot(domain.x, manual_D, label="manual D")
        # display(plot!(domain.x, Dopt_mod, label="Dopt mod"))
        
        plot(domain.x, Dopt_mod, label="Dopt mod")
        plot!(domain.x, solEnt.Dopt, label="Dopt")
        display(title!("Num diffusion at t="*string(round(T_reached, sigdigits=3))))
    else
        u0_temp = solEnt.u_approx[end]
        println("Numerical diffusion negative everywhere")
    end
    println("Going to next time evolution")
end

up_modified = u0_temp
# Final resolution
domain.Tf = T_reached
solEnt = optimize_for_entropy(u0, domain, equation, method)
solRus = fv_solve(domain, u0, equation, Rusanov(CFL_factor))
u_exact = [uexact_burgers_article(xi, solEnt.domain.Tf) for xi in domain.x]

plot(domain.x, u_exact, label="exact")
plot!(domain.x, up_modified, label="modified "*get_name(method))
plot!(domain.x, solEnt.u_approx[end], label=get_name(method))
plot!(domain.x, solRus.u_approx[end], label="Rusanov")