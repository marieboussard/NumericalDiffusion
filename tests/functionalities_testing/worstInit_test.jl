include("../../src/include_file.jl")

# # 1 # With Burgers

# #xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
# CFL_factor = 0.5
# #omega = createInterval(xmin, xmax, Nx, t0, Tf)

# Nx = 100
# equation = burgers()
# #scheme = FVScheme(Euler(), Rusanov(CFL_factor))
# scheme = FVScheme(Euler(), Roe(CFL_factor))


# sol = iterate_WID(-2, 2, Nx, equation, scheme, nb_it=2, boxBounds=[-3 3;])
# #@time sol = iterate_WID(Nx, burgers(), Roe(CFL_factor), nb_it=10, boxBounds=[-3 3;])
# plotWorstWD(sol)

# 2 # Saint-Venant
xmin, xmax = 0, 1
Nx = 5
CFL_factor = 0.5
topoHeight = 0.0
eq = SaintVenant(bump_zb(height=topoHeight), 1e-10)
#scheme = FVScheme(Euler(), Rusanov(CFL_factor))
scheme = FVScheme(Euler(), HR(CFL_factor, Rusanov(CFL_factor)))

boxBounds=[0.1 10;-5.0 5.0]
sourceBounds=[-5.0, 5.0]

sol = iterate_WID(xmin, xmax, Nx, eq, scheme; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)

plotWorstWD(sol, eq)

# # # Reconstruction of the initial data from the optimization results
# # u_init, z = extendInitialDataToLinear(sol, Nx, boxBounds=boxBounds, sourceBounds=sourceBounds)
# u_init, z = extendInitialDataToK(sol, Nx)
# #plot(u_init[:,1], label="water height")
# plot(u_init[:,1] .+ z, label="surface")
# plot!(z, label="z")

# domain = createInterval(xmin, xmax, Nx, 0.0, 0.1)
# @show dt = sol.method.CFL_factor * domain.dx / CFL_cond(sol.equation, u_init) # Timestep given by CFL condition
# # dt = 0.011894258322405449
# @show domain.dx

# domain = createUnitInterval(Nx, 0.0, dt)# Redefining the domain with dt as final time
# domain.sourceVec = z

# fv_sol = fv_solve(domain, u_init, eq, method)
# display(plot_fv_sol(fv_sol, eq, nb_plots=2))

# solEnt = optimize_for_entropy(u_init, domain, eq, method)#, modifiedDataType=maxK())
# display(plot_solution(solEnt))