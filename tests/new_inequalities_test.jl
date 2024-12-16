include("../src/include_file.jl")

xmin, xmax, Nx, t0 = -2, 2, 100, 0
CFL_factor = 0.1
equation = burgers()
scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
#scheme = FVScheme(RK2(), Rusanov(0.5))
#scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod(), domain))

testcase = ArticleTestcase()
#testcase = SimpleShock()
#domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)
Tf = 0.4
domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)
#modifiedDataType = minK()
modifiedDataType = AsymmetricModifiedData()

solEnt = optimize_for_entropy(u0, domain, equation, scheme; modifiedDataType=modifiedDataType)

plot_solution(solEnt)
e1, e2 = checkInequalities(solEnt)

plot(domain.x, e1, label="e1")
display(plot!(domain.x, e2, label="e2"))

display(plot(domain.x, solEnt.u_approx[end]))

@show maximum(e1)
@show maximum(e2)
