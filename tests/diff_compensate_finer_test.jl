include("../src/include_file.jl")

# 1 # Solving Burgers equation

# Domain
xmin, xmax, Nx, t0, Tf = -2, 2, 10, 0, 0.1
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)

# Equation
#equation = burgers()
equation = newEq()

# Finite volume method
# method = Roe(CFL_factor)
method = Centered(CFL_factor)
# method = Rusanov(CFL_factor)

# Initial conditions
testcase = ConcaveConvexTestcase()
#testcase = ArticleTestcase()

# Parameters of diffusion compensation
searchWhenAlreadyEntropic = false

# Profiling time taken by execution
const to = TimerOutput();

function iterate_diffusion_compensation(domain::Domain, equation::Equation, method::FVMethod, testcase::Testcase, searchWhenAlreadyEntropic=false)

    Tf = domain.Tf
    T_reached = 0.0
    timer = 0
    domain_temp = domain # To avoid side effects on the function argument domain

    u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)
    u0_temp = u0

    solAlphaG = OptAlphaGSol{eltype(u0)}()

    while T_reached < Tf

        timer += 1
        dx = domain.dx
        dt = method.CFL_factor * dx / CFL_cond(equation, u0_temp)
        domain_temp.Tf = dt # Tf is modified in the temporary domain so that only one timestep is performed
        T_reached += dt
        solEnt = optimize_for_entropy(u0_temp, domain_temp, equation, method)
        solEntRus = optimize_for_entropy(u0_temp, domain_temp, equation, Rusanov(CFL_factor))

        if maximum(solEnt.Dopt) > 1e-8 || searchWhenAlreadyEntropic==true
            println("Modifiying the scheme to compensate positive numerical diffusion")

            Gent = exactG(Rusanov(CFL_factor), equation, u0_temp)

            # plot(domain.interfaces, solEntRus.Gopt, label="Gopt")
            # display(plot!(domain.interfaces, Gent, label="Gent"))
            # Gent = solEntRus.Gopt

            find_optimal_alphaG(solEnt, Rusanov(CFL_factor), Gent)
            solAlphaG = find_optimal_alphaG(solEnt, Rusanov(CFL_factor), Gent);
            # solAlphaG = find_optimal_alphaG_vec_cons(solEnt, Rusanov(CFL_factor), Gent);

            alpha_opt, Gopt_mod = solAlphaG.alpha, solAlphaG.G

            u0_temp = solAlphaG.up

            plot(domain.interfaces, Gopt_mod, label="Gopt mod")
            plot!(domain.interfaces, solEnt.Gopt, label="Gopt"*get_name(method))
            plot!(domain.interfaces, solEntRus.Gopt, label="Gopt Rusanov")
            display(title!("Num entropy flux at t="*string(round(T_reached, sigdigits=3))))

            Dopt_mod = diffusion(solEnt.u_approx[end-1], u0_temp, Gopt_mod, domain.dx, solEnt.dt_vec[end], equation, domain)
            
            plot(domain.x, Dopt_mod, label="Dopt mod")
            plot!(domain.x, solEnt.Dopt, label="Dopt")
            display(title!("Num diffusion at t="*string(round(T_reached, sigdigits=3))))

            display(plot(domain.interfaces, alpha_opt, label = "alpha opt"))

            alphaG_opt = repasteAlphaG(alpha_opt, Gopt_mod)
            #@show constraint_alphaG(alphaG_opt, solEnt, Rusanov(CFL_factor))
        else
            u0_temp = solEnt.u_approx[end]
            println("Numerical diffusion negative everywhere")
        end
        println("Going to next time evolution")
    end

    println("Total number of iterations: "*string(timer))

    solAlphaG, T_reached

end

function plot_final_results_diffusion_compensation(domain::Domain, equation::Equation, method::FVMethod, testcase::Testcase, solAlphaG::OptAlphaGSol, T_reached)

    alpha_opt = solAlphaG.alpha
    u_modified = solAlphaG.up
    u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)

    # Final resolution
    domain_reached = domain
    domain_reached.Tf = T_reached
    solEnt = optimize_for_entropy(u0, domain_reached, equation, method)
    solRus = fv_solve(domain_reached, u0, equation, Rusanov(CFL_factor))
    u_exact = [uexact_fun(testcase, xi, solEnt.domain.Tf) for xi in domain_reached.x]


    plot(domain_reached.x, u_exact, label="exact")
    plot!(domain_reached.x, u_modified, label="modified "*get_name(method))
    plot!(domain_reached.x, solEnt.u_approx[end], label=get_name(method))
    display(plot!(domain_reached.x, solRus.u_approx[end], label="Rusanov"))

    display(plot(domain_reached.interfaces, alpha_opt, label = "alpha opt"))

end

#@code_warntype iterate_diffusion_compensation(domain, equation, method, testcase, searchWhenAlreadyEntropic)
solAlphaG, T_reached = iterate_diffusion_compensation(domain, equation, method, testcase, searchWhenAlreadyEntropic)
#plot_final_results_diffusion_compensation(domain, equation, method, testcase, solAlphaG, T_reached)
print_timer(to)