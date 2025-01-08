include("../../src/include_file.jl");

Nx = 50
xmin, xmax = -2.0, 2.0
CFL_factor = 0.5
equation = burgers()
testcase = ArticleTestcase()
t0, Tf = 0.0, 0.2
domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)
#domain, u0 = createOneTimestepInterval(Nx, 0.0, xmin, xmax, equation, testcase, CFL_factor)
modifiedDataType = AsymmetricModifiedData()
#scheme = FVScheme(Euler(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), Roe(CFL_factor))
#scheme = FVScheme(RK2(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))

fv_sol = fv_solve(domain, u0, equation, scheme)
plot_fv_sol(fv_sol)
display(title!(get_name(scheme)))

solEntArticle = optimize_for_entropy(u0, domain, equation, scheme; modifiedDataType=modifiedDataType)
@show solEntArticle.optimResult

solEntSimple = optimize_for_entropy(u0, domain, equation, scheme; modifiedDataType=modifiedDataType, boundsType=SimpleBounds())
@show solEntSimple.optimResult

plot(domain.x, solEntArticle.Dopt, label="Article")
plot!(domain.x, solEntSimple.Dopt, label="New Bounds")
display(title!(get_name(scheme)))