include("../src/include_file.jl")

# # 1 # Burgers equation

# Nx = 10

# CFL_number = 0.5
# domain = createInterval(-2, 2, Nx, 0, 0.4)
# u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
# plot(u0)

# @time sol = fv_solve(domain, u0, burgers(), Rusanov(CFL_number))
# plot_fv_sol(sol, uexact_burgers_article)

# @time sol = optimize_for_entropy(u0, domain, burgers(), Rusanov(CFL_number), modifiedDataType=maxK())
# display(plot_solution(sol))
# plot_bounds(sol)


#compare_exact_flux(sol)

# 2 # Saint Venant

Nx, t0, Tf = 50, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
height_bump = 1.0
source_term = bump_zb(height_bump)
eq = SaintVenant(source_term, 1e-10)
addSource!(eq.source, domain)

#eq = SaintVenant(NullSource())

#method = Rusanov(CFL_factor)
method = createHydrostatic(CFL_factor, Rusanov)

v0 = v0_lake_at_rest(domain.x, source_term)
#v0 = v0_lake_at_rest_perturbated(domain.x, source_term)
#v0 = v0_lake_at_rest(domain.x, NullSource())

plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x), label="Water surface")
plot!(domain.x, zb(eq.source, domain.x), label="Topography")

solSV = fv_solve(domain, v0, eq, method)
nb_plots = 5
display(plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots))

solEnt = optimize_for_entropy(v0, domain, eq, method, modifiedDataType=meanK(1,1))
display(plot_solution(solEnt))

plot_bounds(solEnt)