using BenchmarkTools
using NumericalDiffusion

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle
#timescheme = Euler()
timescheme = RK2()

sol = solve(equation, params, timescheme, Rusanov(); log_config=LogConfig(true, false, true, false, false));

estimate = quantify_diffusion(sol, Posteriori(AsymmetricMD()));

estimator = Estimator(sol, Posteriori(AsymmetricMD()), 0);
eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
eta!(estimator.entfun, estimator.u, estimator.etacont)
compute_G_bounds!(estimator)

etavec = zero(mesh.x)
Gvec = zero(mesh.x)
eta!(BurgersEnt(), estimate.uinit, etavec)
G!(BurgersEnt(), estimate.uinit, Gvec)
Gexact = zero(mesh.x)
for i in 1:Nx
    Gexact[i] = 0.5*(Gvec[i] + Gvec[mod1(i+1,Nx)]) - max(abs(estimate.uinit[i]), abs(estimate.uinit[mod1(i+1,Nx)]))*0.5*(etavec[mod1(i+1,Nx)] - etavec[i])
end

using Plots
plot(mesh.x, sol.uinit, label="uinit")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))
plot(mesh.x, estimate.m, label="m")
plot!(mesh.x, estimate.M, label="M")
#plot!(mesh.x, Gexact, label="Gexact")
display(plot!(mesh.x, estimate.Gopt, label="Optimal Numerical Entropy Flux"))
display(plot(mesh.x, estimate.D, label="Numerical Diffusion"))
# axislegend(position = :rb)
# current_figure()

#@btime estimate = quantify_diffusion(sol, Posteriori(AsymmetricMD()));