# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 1.1 # Solving Burgers equation with Rusanov

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0.0, 0.4
CFL_factor = 0.5
domain = Interval(Nx, xmin, xmax, t0, Tf)
equation = burgers()
testcase = ArticleTestcase()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
saveLog = true

u_init = initialData(domain, testcase)
uexact = exactData(domain, testcase)

problem = FVProblem(domain, equation, testcase, scheme, saveLog, u_init)
@time fv_sol = solve(problem)
plot_fv_sol(fv_sol, 4)
display(title!("Rusanov"))

# solBurgers = fv_solve(domain, u0, burgers(), scheme)
# #plot_fv_sol(solBurgers, nb_plots=6)
# plot_fv_sol(solBurgers, uexact_burgers_article)
# display(title!("Rusanov"))


# 1.2 # Solving Burgers equation with Roe

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0.0, 0.4
CFL_factor = 0.5
domain = Interval(Nx, xmin, xmax, t0, Tf)
equation = burgers()
testcase = ArticleTestcase()
scheme = FVScheme(Euler(), Roe(CFL_factor))
saveLog = true

u_init = initialData(domain, testcase)
uexact = exactData(domain, testcase)

problem = FVProblem(domain, equation, testcase, scheme, saveLog, u_init)
fv_sol = solve(problem)
plot_fv_sol(fv_sol, 4)
display(title!("Roe"))

# solBurgers = fv_solve(domain, u0, burgers(), scheme)
# plot_fv_sol(solBurgers, nb_plots=6)
# #plot_fv_sol(solBurgers, uexact_burgers_article)
# title!("Roe")