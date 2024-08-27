include("../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 10, 0, 0.4
CFL_factor = 0.5
domain = createInterval(xmin, xmax, Nx, t0, Tf)
u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
equation = burgers()
method = Roe(CFL_factor)

# solBurgers = fv_solve(domain, u0, equation, method)
# plot_fv_sol(solBurgers, uexact_burgers_article)

solEnt = optimize_for_entropy(u0, domain, equation, method)
plot_solution(solEnt)

# g_delta = B_delta(solEnt)
# display(plot(domain.interfaces, g_delta))

# D_test = diffusion(solEnt.u_approx[end-1], solEnt.u_approx[end], solEnt.Gopt .+ g_delta, solEnt.domain.dx, solEnt.dt_vec[end], equation, domain)
# plot(domain.x, D_test)
# plot!(domain.x, solEnt.Dopt./10)

solA = find_optimal_A(solEnt);