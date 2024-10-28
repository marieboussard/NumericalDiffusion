include("../src/include_file.jl")

# 1 # Solving Burgers equation
xmin, xmax, Nx, t0, Tf = -2, 2, 20, 0, 0.4

CFL_factor = 0.5
domain = createInterval(xmin, xmax, Nx, t0, Tf)
#equation = burgers()
equation = newEq()
# method = Roe(CFL_factor)
method = Centered(CFL_factor)
#method = Rusanov(CFL_factor)
#u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
#u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_new(domain.x[i])] end; res)

testcase = ConcaveConvexTestcase()
#testcase = ArticleTestcase()
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)

T_reached = 0.0
searchWhenAlreadyEntropic = true

u0_temp = u0
for k in 1:3
    global u0_temp
    global T_reached
    global solEnt
    global alpha_opt
    global solEntRus

    dx = domain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u0_temp)
    domain.Tf = dt # Tf is modified so that only one timestep is performed
    T_reached += dt
    solEnt = optimize_for_entropy(u0_temp, domain, equation, method)
    solEntRus = optimize_for_entropy(u0_temp, domain, equation, Rusanov(CFL_factor))
    # plot_solution(solEntRus)

    if maximum(solEnt.Dopt) > 1e-8 || searchWhenAlreadyEntropic==true
        println("Modifiying the scheme to compensate positive numerical diffusion")

        # Gent = exactG(Rusanov(CFL_factor), equation, u0_temp)

        # plot(domain.interfaces, solEntRus.Gopt, label="Gopt")
        # display(plot!(domain.interfaces, Gent, label="Gent"))
        Gent = solEntRus.Gopt

        # solAlphaG = find_optimal_alphaG(solEnt, Rusanov(CFL_factor), Gent);
        solAlphaG = find_optimal_alphaG_vec_cons(solEnt, Rusanov(CFL_factor), Gent);

        alpha_opt, Gopt_mod = solAlphaG.alpha, solAlphaG.G

        u0_temp = solAlphaG.up
        #u0_temp = up_mod

        #u0_temp = modified_scheme_step(Aopt, solEnt)
        # plot(domain.x, u0_temp, label="u0_temp")
        # display(plot!(domain.x, up_mod, label="up mod"))


        # plot(domain.interfaces, m_delta, label="m delta")
        # plot!(domain.interfaces, Gopt_mod, label="Gopt mod")
        # plot!(domain.interfaces, M_delta, label="M delta")
        # display(title!("Bounds at t="*string(round(T_reached, sigdigits=3))))

        plot(domain.interfaces, Gopt_mod, label="Gopt mod")
        plot!(domain.interfaces, solEnt.Gopt, label="Gopt"*get_name(method))
        plot!(domain.interfaces, solEntRus.Gopt, label="Gopt Rusanov")
        display(title!("Num entropy flux at t="*string(round(T_reached, sigdigits=3))))

        Dopt_mod = diffusion(solEnt.u_approx[end-1], u0_temp, Gopt_mod, domain.dx, solEnt.dt_vec[end], equation, domain)
        
        # plot(domain.x, manual_D, label="manual D")
        # display(plot!(domain.x, Dopt_mod, label="Dopt mod"))
        
        plot(domain.x, Dopt_mod, label="Dopt mod")
        plot!(domain.x, solEnt.Dopt, label="Dopt")
        display(title!("Num diffusion at t="*string(round(T_reached, sigdigits=3))))

        display(plot(domain.interfaces, alpha_opt, label = "alpha opt"))

        alphaG_opt = repasteAlphaG(alpha_opt, Gopt_mod)
        @show constraint_alphaG(alphaG_opt, solEnt, Rusanov(CFL_factor))
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
# u_exact = [uexact_burgers_article(xi, solEnt.domain.Tf) for xi in domain.x]
# u_exact = [uexact_newEq(xi, solEnt.domain.Tf) for xi in domain.x]
u_exact = [uexact_fun(testcase, xi, solEnt.domain.Tf) for xi in domain.x]


plot(domain.x, u_exact, label="exact")
plot!(domain.x, up_modified, label="modified "*get_name(method))
plot!(domain.x, solEnt.u_approx[end], label=get_name(method))
display(plot!(domain.x, solRus.u_approx[end], label="Rusanov"))

plot(domain.interfaces, alpha_opt, label = "alpha opt")
