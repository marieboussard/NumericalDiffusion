include("../../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
testcase = ArticleTestcase()
u0 = initialData(domain, testcase)
equation = burgers()
scheme = FVScheme(Euler(), Rusanov(CFL_factor))

solEnt = optimize_for_entropy(u0, domain, equation, scheme)
solEnt_light = optimize_for_entropy(u0, domain, equation, scheme; boundsType=LightBounds())

plot(domain.x, solEnt.Dopt, label="Normal Bounds")
plot!(domain.x, solEnt_light.Dopt, label="Light Bounds")

# # 1 # Hydrostatic scheme for the lake at rest

# Nx, t0, Tf = 100, 0, 0.4
# CFL_factor = 0.5
# height = 0.5
# width = 0.5
# topography = bump_zb(height=height, width=width)
# domain = createUnitInterval(Nx, t0, Tf)
# eq = SaintVenant(topography, 1e-10)
# addSource!(eq.source, domain)
# method = createHydrostatic(CFL_factor, Rusanov)
# #method = Rusanov(CFL_factor)

# v0 = v0_lake_at_rest(domain.x, topography)

# solEnt = optimize_for_entropy(v0, domain, eq, method; boundsType=LightBounds())

# plot_solution(solEnt)