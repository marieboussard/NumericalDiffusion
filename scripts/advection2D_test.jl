using FiniteVolumes
using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
Ny = 100
xmin, xmax, ymin, ymax = -2, 2,-2,2
t0, tf = 0.0, 1.4
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

#equation = Advection2Example
equation = advection2_vecfield(mesh)

sol = solve(equation, params, Euler(), Rusanov2D(); log_config=LogConfig(true,true,true,true));

clim = (minimum(sol.uinit), maximum(sol.uinit)) # Fixing color range
display(heatmap(mesh.x, mesh.y, sol.uinit, clim=clim, aspect_ratio=:equal, title="initial"))
heatmap(mesh.x, mesh.y, sol.u, clim=clim, aspect_ratio=:equal, title="final")