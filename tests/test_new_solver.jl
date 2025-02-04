include("../src/include_file.jl")

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0.0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = ArticleTestcase()
u_init = initialData(domain, testcase)
uexact = exactData(domain, testcase)
equation = burgers()
saveLog = true

scheme = FVScheme(Euler(), Rusanov(CFL_factor))

problem = FVProblem(domain, equation, testcase, scheme, saveLog, u_init)

plot(domain.x, u_init, label="init")

fv_sol = solve(problem)
plot_fv_sol(fv_sol, 4)