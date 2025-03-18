using BenchmarkTools
include("../../src/numdiff/include_file.jl")

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false));

estimate_priori = quantify_diffusion(sol, PrioriMultidim(); name="priori multidim");

using Plots
plot(mesh.x, estimate.uinit, label="uinit")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))
plot(mesh.x, estimate.l, label="l")
display(plot!(mesh.x, estimate.L, label="L"))
plot(mesh.x, estimate.Dlow, label="Dlow")
plot(mesh.x, estimate.D, label="Numerical Diffusion")