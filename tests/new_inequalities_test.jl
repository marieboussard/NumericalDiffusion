include("../src/include_file.jl")

xmin, xmax, Nx, t0 = -2, 2, 20, 0
CFL_factor = 0.5
equation = burgers()
#method = Roe(CFL_factor)
method = Rusanov(CFL_factor)
testcase = ArticleTestcase()
#testcase = SimpleShock()
domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)
modifiedDataType = minK()

solEnt = optimize_for_entropy(u0, domain, equation, method; modifiedDataType=modifiedDataType)

plot_solution(solEnt)
e1, e2 = checkInequalities(solEnt)

plot(domain.x, e1, label="e1")
plot!(domain.x, e2, label="e2")
