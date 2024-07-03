#using Plots

#include("../src/tools/fv_solution.jl")
#include("../src/tools/method.jl")
#include("../src/tools/equation.jl")
#include("../src/tools/domain.jl")
#include("../src/opt_diffusion.jl")

include("../src/include_file.jl")

# # 1 # Burgers equation

# Nx = 100

# CFL_number = 0.5
# domain = createInterval(-2, 2, Nx, 0, 0.4)
# u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(domain.x[i])] end; res)
# plot(u0)

# @time sol = fv_solve(domain, u0, burgers(), Rusanov(CFL_number))
# plot_fv_sol(sol, uexact_burgers_article)

# @time sol = optimize_for_entropy(u0, domain, burgers(), Rusanov(CFL_number), modifiedDataType=maxK())
# plot_solution(sol)
# #plot_bounds(sol)


#compare_exact_flux(sol)

# 2 # Saint Venant

Nx, t0, Tf = 5, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
eq = SaintVenant(Bump_zb())
#eq = SaintVenant(NullSource())

v0 = v0_lake_at_rest(domain.x, Bump_zb())
#v0 = v0_lake_at_rest(domain.x, NullSource())

plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x), label="Water surface")
plot!(domain.x, zb(eq.source, domain.x), label="Topography")

# solSV = fv_solve(domain, v0, eq, Rusanov(CFL_factor))

# nb_plots = 5

# plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)

solEnt = optimize_for_entropy(v0, domain, eq, Rusanov(CFL_number))
plot_solution(solEnt)