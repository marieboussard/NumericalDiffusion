include("../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
CFL_factor = 0.5
domain = createInterval(xmin, xmax, Nx, t0, Tf)
#global u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
equation = burgers()
method = Roe(CFL_factor)
#method = Rusanov(CFL_factor)
# dx = domain.dx
# dt = method.CFL_factor * dx / CFL_cond(equation, u0)
# domain.Tf = dt # Tf is modified so that only one timestep is performed
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)

T_reached = 0.0
gmidShape = ConsistencyMeanShape()
#gmidShape = RusanovLikeShape()
searchWhenAlreadyEntropic = true

# solBurgers = fv_solve(domain, u0, equation, method)
# plot_fv_sol(solBurgers, uexact_burgers_article)
u0_temp = u0
for k in 1:1
    global u0_temp
    global T_reached
    global solEnt
    global Aopt
    global solEntRus

    dx = domain.dx
    dt = method.CFL_factor * dx / CFL_cond(equation, u0_temp)
    domain.Tf = dt # Tf is modified so that only one timestep is performed
    T_reached += dt
    solEnt = optimize_for_entropy(u0_temp, domain, equation, method)
    #solRus = fv_solve(domain, u0, equation, Rusanov(CFL_factor))
    solEntRus = optimize_for_entropy(u0_temp, domain, equation, Rusanov(CFL_factor))

    #=
    @show A_init, _p = zero(solEnt.Gopt), 0.0
    cons(res, A, p) = (res .= [constraint_A(A, sol; gmidShape=gmidShape)])
    @show constraint_A(A_init, solEnt; gmidShape=gmidShape)

end
=#


    if maximum(solEnt.Dopt) > 1e-8 || searchWhenAlreadyEntropic==true
        println("Modifiying the scheme to compensate positive numerical diffusion")

        @enter solA = find_optimal_A(solEnt; gmidShape=gmidShape)
        Aopt = solA.u
        u0_temp = modified_scheme_step(Aopt, solEnt)

        Gmid = compute_G_mid(Aopt, solEnt, gmidShape)
        plot(domain.interfaces, Gmid, label="G mid")
        plot!(domain.interfaces, solEnt.Gopt, label="Gopt"*get_name(method))
        plot!(domain.interfaces, solEntRus.Gopt, label="Gopt Rusanov")
        display(title!("t="*string(round(T_reached, sigdigits=3))))

        Dmid = diffusion_from_A(Aopt, solEnt; gmidShape=gmidShape)
        plot(domain.x, Dmid, label="Dmid")
        plot!(domain.x, solEnt.Dopt, label="Dopt")
        display(title!("t="*string(round(T_reached, sigdigits=3))))
    else
        u0_temp = solEnt.u_approx[end]
        println("Numerical diffusion negative everywhere")
    end
    println("Going to next time evolution")
end

@show T_reached
#plot_solution(solEnt)

# g_delta = B_delta(solEnt)
# display(plot(domain.interfaces, g_delta))

# D_test = diffusion(solEnt.u_approx[end-1], solEnt.u_approx[end], solEnt.Gopt .+ g_delta, solEnt.domain.dx, solEnt.dt_vec[end], equation, domain)
# plot(domain.x, D_test)
# plot!(domain.x, solEnt.Dopt./10)

#solA = find_optimal_A(solEnt);
#Aopt = solA.u
#up_modified = modified_scheme_step(Aopt, u0_temp)
up_modified = u0_temp
# Final resolution
domain.Tf = T_reached
solEnt = optimize_for_entropy(u0, domain, equation, method)
solRus = fv_solve(domain, u0, equation, Rusanov(CFL_factor))
u_exact = [uexact_burgers_article(xi, solEnt.domain.Tf) for xi in domain.x]

# Gmid = compute_G_mid(Aopt, solEnt)
# plot(domain.interfaces, Gmid, label="G mid")
# display(plot!(domain.interfaces, solEnt.Gopt, label="Gopt"))

plot(domain.x, u_exact, label="exact")
plot!(domain.x, up_modified, label="modified "*get_name(method))
plot!(domain.x, solEnt.u_approx[end], label=get_name(method))
plot!(domain.x, solRus.u_approx[end], label="Rusanov")


