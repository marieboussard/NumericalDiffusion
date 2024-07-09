include("../src/include_file.jl")

Nx, t0, Tf = 10, 0.0, 0.1
sourceHeight = 1.0
CFL_factor = 0.5

domain = createUnitInterval(Nx, t0, Tf)
eq = SaintVenant(Bump_zb(sourceHeight), 1e-10)
#method = createHydrostatic(CFL_factor, Rusanov)
method = Rusanov(CFL_factor)
addSource!(eq.source, domain)

v0 = v0_lake_at_rest(domain.x, eq.source)

# plot(domain.x, v0[:,1] .+ domain.sourceVec)
# plot!(domain.x, domain.sourceVec)

solSV = fv_solve(domain, v0, eq, method)

nb_plots = 5
plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)