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

estim_prio = quantify_diffusion(sol, Priori(AsymmetricMD()));
estim_priomulti = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));

using Plots
plot(mesh.x, estim_prio.D, label="priori")
plot!(mesh.x, estim_priomulti.D, label = "priori mulitidim")