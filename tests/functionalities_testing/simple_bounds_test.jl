include("../../src/include_file.jl")

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
testcase = ArticleTestcase()

domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)
#scheme = FVScheme(Euler(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), Roe(CFL_factor))
#scheme = FVScheme(RK2(), Rusanov(CFL_factor))
scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
equation = burgers()

@time solSimple = optimize_for_entropy(u0, domain, equation, scheme, boundsType=SimpleBounds())
@time sol = optimize_for_entropy(u0, domain, equation, scheme, modifiedDataType=AsymmetricModifiedData())

plot_solution(solSimple)
plot_solution(sol)