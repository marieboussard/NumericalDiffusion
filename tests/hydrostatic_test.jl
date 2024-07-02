include("../src/include_file.jl")

# # 1 # Hydrostatic scheme for the lake at rest

# Nx, t0, Tf = 100, 0, 0.4
# CFL_factor = 0.5
# domain = createUnitInterval(Nx, t0, Tf)
# eq = SaintVenant(Bump_zb())

# v0 = v0_lake_at_rest(domain.x, Bump_zb())

# plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
# plot!(domain.x, zb(eq.source, domain.x))

# solSV = fv_solve(domain, v0, eq, createHydrostatic(CFL_factor, Rusanov))

# nb_plots = 5

# plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)


# 2 # Hydrostatic scheme for the perturbated lake at rest

Nx, t0, Tf = 100, 0, 0.1
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
eq = SaintVenant(Bump_zb())

v0 = v0_lake_at_rest_perturbated(domain.x, Bump_zb())

plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
plot!(domain.x, zb(eq.source, domain.x))

solSV = fv_solve(domain, v0, eq, createHydrostatic(CFL_factor, Rusanov))

nb_plots = 5

plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)