include("../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)
eq = burgers()
method = Rusanov(CFL_factor)

solEnt = optimize_for_entropy(u0, domain, eq, method; boundsType=LightBounds())
plot_solution(solEnt)

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