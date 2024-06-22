using Plots

include("../src/tools/fv_solution.jl")
include("../src/tools/method.jl")
include("../src/tools/equation.jl")
include("../src/tools/domain.jl")

CFL_number = 0.5
domain = createInterval(-2, 2, 100, 0, 0.4)
sol = fv_solve(domain, u0_burgers_article, Burgers(), Rusanov(CFL_number))
plot_fv_sol(sol, uexact_burgers_article)