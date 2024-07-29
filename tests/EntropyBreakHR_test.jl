include("../src/include_file.jl")

Nx = 50
xmin, xmax = 0.0, 1.0
CFL_factor = 0.5
source = Discontinuous_zb()
eq = SaintVenant(source, 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
domain = createInterval(xmin, xmax, Nx, 0.0, 0.1)
dt = sol.method.CFL_factor * domain.dx / CFL_cond(sol.equation, u_init)
domain = createInterval(xmin, xmax, Nx, 0.0, dt)
addSource!(source, domain)
u_init = v0_discontinuous(domain.x)

fv_sol = fv_solve(domain, u_init, eq, method)
display(plot_fv_sol(fv_sol, eq, nb_plots=2))
solEnt = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, g_tol=1e-10)#, modifiedDataType=maxK())
plot_solution(solEnt)