include("../src/include_file.jl")

# Solving Burgers equation with MUSCL (constructed with Rusanov)

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = ArticleTestcase()
u0 = initialData(domain, testcase)
#scheme = FVScheme(Euler(), MUSCL(0.5, Rusanov(0.5), Minmod(), domain))
#scheme = FVScheme(RungeKutta(), Rusanov(0.5))
scheme = FVScheme(RungeKutta(), MUSCL(0.1, Rusanov(0.1), Minmod(), domain))

solBurgers_muscl = fv_solve(domain, u0, burgers(), scheme)
solBurgers_Rusanov = fv_solve(domain, u0, burgers(), FVScheme(Euler(), Rusanov(0.5)))
#plot_fv_sol(solBurgers, nb_plots=6)
plot_fv_sol(solBurgers_muscl, uexact_burgers_article)
plot!(domain.x, solBurgers_Rusanov.u_approx[end], label="Rusanov")
