include("../src/include_file.jl")

# 1.1 # Solving "newEq" equation with Rusanov

xmin, xmax, Nx, t0, Tf = -2, 2, 60, 0, 0.3
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)
testcase = ConcaveConvexTestcase()
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, omega.x[i])] end; res)
#u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_new(omega.x[i])] end; res)


solNew = fv_solve(omega, u0, newEq(), Rusanov(CFL_factor))
#plot_fv_sol(solBurgers, nb_plots=6)

u_exact = [uexact_fun(testcase, xi, omega.Tf) for xi in omega.x]
#plot_fv_sol(solNew, nb_plots=5)
plot_fv_sol(solNew, (x,t) -> uexact_fun(testcase, x, t))
display(title!("Rusanov"))

# # Quantifying numerical diffusion

# sol = optimize_for_entropy(u0, omega, newEq(), Rusanov(CFL_factor))

# display(plot_solution(sol))
# @show sol.Jopt