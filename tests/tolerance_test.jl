include("../src/include_file.jl")

Nx = 50
CFL_factor = 0.5
eq = SaintVenant(bump_zb(), 1e-10)

method = createHydrostatic(CFL_factor, Rusanov)

# boxBounds=[0.0 3;-1.0 1.0]
# sourceBounds=[-1.0, 1.0]

# # addSource!(eq.source, domain)
# # u_init = v0_lake_at_rest(domain.x, eq.source)

# sol = iterate_WID(Nx, eq, method; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)

# # Reconstruction of the initial data from the optimization results
# u_init, z = extendInitialDataToK(sol, Nx)

# domain = createUnitInterval(Nx, 0.0, 0.1)
# dt = sol.method.CFL_factor * domain.dx / CFL_cond(sol.equation, u_init) # Timestep given by CFL condition

# # Redefining the domain with dt as final time
# domain = createUnitInterval(Nx, 0.0, dt)
# domain.sourceVec = z

# # Solving Saint-Venant for this data
# fv_sol = fv_solve(domain, u_init, eq, method)
# display(plot_fv_sol(fv_sol, eq, nb_plots=2))

# Optimization process

solEnt = optimize_for_entropy(u_init, domain, eq, method; f_tol=1e-12, iterations=100000)#, g_tol=1e-10)#, modifiedDataType=maxK())
plot_solution(solEnt)