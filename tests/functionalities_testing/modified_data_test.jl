""" Testing different kinds of modified data on Saint Venant equations """

include("../src/include_file.jl")

# 1 # Hydrostatic scheme for the lake at rest

Nx, t0, Tf = 100, 0, 0.4
CFL_factor = 0.5
height = 0.5
width = 0.5
topography = bump_zb(height=height, width=width)
domain = createUnitInterval(Nx, t0, Tf)
eq = SaintVenant(topography, 1e-10)
addSource!(eq.source, domain)
method = createHydrostatic(CFL_factor, Rusanov)
#method = Rusanov(CFL_factor)

v0 = v0_lake_at_rest(domain.x, topography)

# plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
# plot!(domain.x, zb(eq.source, domain.x))

# solSV = fv_solve(domain, v0, eq, createHydrostatic(CFL_factor, Rusanov))

# nb_plots = 5
# plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)
solEntExtend = optimize_for_entropy(v0, domain, eq, method; modifiedDataType=AsymmetricModifiedData())
solEntK = optimize_for_entropy(v0, domain, eq, method)
solEntExtend.label = "extend"
solEntK.label="K"

display(plot_solutions([solEntExtend, solEntK]))