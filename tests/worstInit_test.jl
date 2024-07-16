include("../src/include_file.jl")

# # 1 # With Burgers

# #xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
# CFL_factor = 0.5
# #omega = createInterval(xmin, xmax, Nx, t0, Tf)

# Nx = 100

# #@time sol = iterate_WID(Nx, burgers(), Rusanov(CFL_factor), nb_it=10, boxBounds=[-3 3;])
# @time sol = iterate_WID(Nx, burgers(), Roe(CFL_factor), nb_it=10, boxBounds=[-3 3;])
# plotWorstWD(sol)

# 2 # Saint-Venant
Nx = 50
CFL_factor = 0.5
topoHeight = 2.0
eq = SaintVenant(Bump_zb(topoHeight), 1e-10)

method = createHydrostatic(CFL_factor, Rusanov)
#method = Rusanov(CFL_factor)

boxBounds=[0.0 3;-1.0 1.0]
sourceBounds=[-1.0, 1.0]

sol = iterate_WID(Nx, eq, method; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)

plotWorstWD(sol, eq)

# # Reconstruction of the initial data from the optimization results
u_init, z = extendInitialDataToLinear(sol, Nx, boxBounds=boxBounds, sourceBounds=sourceBounds)
# # plot(u_init[:,1], label="water height")
# # plot!(u_init[:,1] .+ z, label="surface")
# # plot!(z, label="z")

domain = createUnitInterval(Nx, 0.0, 0.1)
domain.sourceVec = z

fv_sol = fv_solve(domain, u_init, eq, method)
display(plot_fv_sol(fv_sol, eq, nb_plots=5))

solEnt = optimize_for_entropy(u_init, domain, eq, method)#, modifiedDataType=maxK())
display(plot_solution(solEnt))