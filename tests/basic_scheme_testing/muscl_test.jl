include("../../src/include_file.jl")

println("===========BEGINNING===============")

# Solving Burgers equation with MUSCL (constructed with Rusanov)

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
testcase = ArticleTestcase()

domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)

#domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, burgers(), testcase, CFL_factor)

#scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
#scheme = FVScheme(RK2(), Rusanov(0.5))
# scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Superbee()))

# 1 ## Applying the scheme

println("MUSCL")
solBurgers_muscl = fv_solve(domain, u0, burgers(), scheme)
println("Rusanov")
solBurgers_Rusanov = fv_solve(domain, u0, burgers(), FVScheme(Euler(), Rusanov(CFL_factor)))
#plot_fv_sol(solBurgers, nb_plots=6)
display(plot_fv_sol(solBurgers_muscl, uexact_burgers_article))
display(plot!(domain.x, solBurgers_Rusanov.u_approx[end], label="Rusanov"))

# plot(domain.x[1:2], solBurgers_muscl.u_approx[end][1:2])
# plot!(domain.x[1:2], solBurgers_Rusanov.u_approx[end][1:2], label="Rusanov")

## 2 ## Quantifying Numerical Diffusion

#modifiedDataType = maxK()
modifiedDataType = AsymmetricModifiedData()

sol = optimize_for_entropy(u0, domain, burgers(), scheme, modifiedDataType=modifiedDataType)
@show sol.Jopt
display(plot_solution(sol))

solRus = optimize_for_entropy(u0, domain, burgers(), FVScheme(Euler(), Rusanov(CFL_factor)), modifiedDataType=modifiedDataType)
@show solRus.Jopt
display(plot_solution(solRus))