include("../src/include_file.jl")

# 1 # Solving Burgers equation

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)

solBurgers = fv_solve(omega, u0, burgers(), Rusanov(CFL_factor))
plot_fv_sol(solBurgers, uexact_burgers_article)
#plot_fv_sol(solBurgers, nb_plots=6)