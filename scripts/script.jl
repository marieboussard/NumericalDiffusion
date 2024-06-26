using Plots

#include("../src/tools/fv_solution.jl")
#include("../src/tools/method.jl")
#include("../src/tools/equation.jl")
#include("../src/tools/domain.jl")
include("../src/opt_diffusion.jl")

CFL_number = 0.5
domain = createInterval(-2, 2, 100, 0, 0.4)
@time sol = fv_solve(domain, u0_burgers_article, burgers(), Rusanov(CFL_number))
plot_fv_sol(sol, uexact_burgers_article)

# @time sol = optimize_for_entropy(u0_burgers_article, domain, Burgers(), Roe(CFL_number))
# plot_solution(sol)
#plot_bounds(sol)


#compare_exact_flux(sol)