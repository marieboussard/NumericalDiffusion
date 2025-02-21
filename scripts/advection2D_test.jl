using FiniteVolumes
using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 10
Ny = 10
xmin, xmax, ymin, ymax = -2, 2,-2,2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

equation = Advection2Example

sol = solve(equation, params, Euler(), Rusanov2D(); log_config=LogConfig(true,true,true,true));

display(heatmap(mesh.x, mesh.y, sol.uinit, title="initial"))
heatmap(mesh.x, mesh.y, sol.u, title="final")