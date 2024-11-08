# A serie of tests that need to be checked at each major modification of the code

include("../../src/include_file.jl")

# 1.1 # Solving Burgers equation with Rusanov

xmin, xmax, Nx, t0, Tf = -2, 2, 15, 0, 0.4
CFL_factor = 0.5
omega = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)

solBurgers = fv_solve(omega, u0, burgers(), Rusanov(CFL_factor))
#plot_fv_sol(solBurgers, nb_plots=6)
plot_fv_sol(solBurgers, uexact_burgers_article)
display(title!("Rusanov"))

@show size(solBurgers.u_approx[end])

# 1.2 # Solving Burgers equation with Roe

xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.4
CFL_factor = 0.5
omega = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)

solBurgers = fv_solve(omega, u0, burgers(), Roe(CFL_factor))
#plot_fv_sol(solBurgers, nb_plots=6)
plot_fv_sol(solBurgers, uexact_burgers_article)
title!("Roe")