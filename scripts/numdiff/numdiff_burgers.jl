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

sol = solve(equation, params, Euler(), Rusanov(); maxiter=1);

estimate = quantify_diffusion(sol, Posteriori());

using Plots
plot(mesh.x, sol.uinit, label="uinit")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))
plot(mesh.x, estimate.m, label="m")
plot!(mesh.x, estimate.M, label="M")
display(plot!(mesh.x, estimate.Gopt, label="Optimal Numerical Entropy Flux"))
plot(mesh.x, estimate.D, label="Numerical Diffusion")