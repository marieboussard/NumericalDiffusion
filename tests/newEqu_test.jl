include("../src/include_file.jl")

# 1.1 # Solving "newEq" equation with Rusanov

xmin, xmax, Nx, t0, Tf = -0.5, 0.5, 15, 0, 0.4
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_new(omega.x[i])] end; res)

solNew = fv_solve(omega, u0, newEq(), Rusanov(CFL_factor))
#plot_fv_sol(solBurgers, nb_plots=6)
plot_fv_sol(solNew, nb_plots=5)
display(title!("Rusanov"))