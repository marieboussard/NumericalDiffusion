include("../src/include_file.jl")

# # 1 # Solving Burgers equation

# xmin, xmax, Nx, t0, Tf = -2, 2, 5, 0, 0.4
# CFL_factor = 0.5
# omega = createInterval(xmin, xmax, Nx, t0, Tf)

# #u0 = u0_burgers_article.(omega.x)
# u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)
# #plot(u0)

# solBurgers = fv_solve(omega, u0, burgers(), Rusanov(CFL_factor))
# plot_fv_sol(solBurgers, uexact_burgers_article)
# #plot_fv_sol(solBurgers, nb_plots=6)


# 2 # Solving Saint Venant equation

Nx, t0, Tf = 20, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
eq = SaintVenant(Bump_zb())

v0 = v0_lake_at_rest(domain.x, Bump_zb())

# plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
# plot!(domain.x, zb(eq.source, domain.x))

solSV = fv_solve(domain, v0, eq, Rusanov(CFL_factor))

nb_plots = 5

plot_fv_sol(solSV, solSV.equation; nb_plots=nb_plots)