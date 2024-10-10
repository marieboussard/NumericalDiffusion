# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 4 # Quantifying numerical diffusion of lake at rest perturbated solutions of Saint-Venant with hydrostatic reconstruction scheme

Nx, t0, Tf = 100, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
topography = bump_zb(height=0.5, width = 0.4)
eq = SaintVenant(topography, 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)

v0 = v0_lake_at_rest_perturbated(domain.x, eq.source)
fv_sol = fv_solve(domain, v0, eq, method)
display(plot_fv_sol(fv_sol, eq; nb_plots=5))
solEnt = optimize_for_entropy(v0, domain, eq, method)
plot_solution(solEnt)