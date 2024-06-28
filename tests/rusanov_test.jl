include("../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)

u0 = u0_burgers_article.(omega.x)
solBurgers = fv_solve(omega, u0, burgers(), Rusanov(CFL_factor))
#plot_fv_sol(solBurgers, uexact_burgers_article)
plot_fv_sol(solBurgers, nb_plots=6)


# # 2 # Solving Saint Venant equation

# Nx, t0, Tf = 10, 0, 0.4
# CFL_factor = 0.5
# domain = createUnitInterval(Nx, t0, Tf)

# v0 = v0_lake_at_rest(domain.x, Bump_zb())

# plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
# plot!(domain.x, zb(eq.source, domain.x))

# solSV = fv_solve(domain, v0, SaintVenant(Bump_zb()), Rusanov(CFL_factor))

# nb_plots = 5
# p = div(sol.Nt, nb_plots)

# plt3 = plot()

# for k in 0:nb_plots-1
#     plot!(domain.x, [solSV.u_approx[k*p+1][i][1] for i in 1:Nx] .+ zb(Bump_zb(), domain.x), label="t = " * string(round(solSV.t_vec[k*p+1], sigdigits=2)))
# end
# plot!(domain.x, zb(Bump_zb(), domain.x), label="zb")
# xlabel!("x")
# display(ylabel!("Surface of the lake"))