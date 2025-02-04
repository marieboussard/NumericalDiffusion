include("../src/include_file.jl")

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = ArticleTestcase()
u0 = initialData(domain, testcase)
uexact = exactData(domain, testcase)
equation = burgers()

scheme = FVScheme(Euler(), Rusanov(CFL_factor))

problem = FVProblem(domain, equation, testcase, scheme, u0)
fv_sol = solve(problem)
plot_fv_sol(fv_sol)

# plot(domain.x, fv_sol.problem.u_init, label="t=0")
# plot!(domain.x, fv_sol.u_approx, label="t="*string(fv_sol.t))
# plot!(domain.x, uexact, label="exact")