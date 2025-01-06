# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 3 # Solving lake at rest states with hydrostatic reconstruction scheme

Nx, t0, Tf = 50, 0, 0.4
CFL_factor = 0.5
height = 1.0
domain = createUnitInterval(Nx, t0, Tf)
#topography = bump_zb(height=0.5, width = 0.4)
#topography = flat_zb()
topography = sinus_zb()
eq = SaintVenant(topography, 1e-10)
addSource!(eq.source, domain)

#scheme = FVScheme(Euler(), Rusanov(CFL_factor))
scheme = FVScheme(Euler(), HR(CFL_factor, Rusanov(CFL_factor)))

v0 = v0_lake_at_rest(domain.x, topography)

plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
plot!(domain.x, zb(eq.source, domain.x))

solSV = fv_solve(domain, v0, eq, scheme)

nb_plots = 5

plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)