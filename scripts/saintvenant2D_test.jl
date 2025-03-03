using FiniteVolumes
using BenchmarkTools
using UnPack

# Domain definition
Nx = 20
Ny = 20
xmin, xmax, ymin, ymax = -2, 2,-2,2
t0, tf = 0.0, 2.0
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

using FiniteVolumes: SaintVenant2Flat
#equation = SaintVenant2Flat
#equation = SaintVenantFlat2
equation = SaintVenantAtRest2

sol = solve(equation, params, Euler(), Rusanov2D(); maxiter=1);

using Plots
clim = (minimum(sol.uinit), maximum(sol.uinit))
display(heatmap(mesh.x, mesh.y, sol.u[:,:,1], clim=clim, aspect_ratio=:equal, title="h"))
display(heatmap(mesh.x, mesh.y, sol.u[:,:,2], clim=clim, aspect_ratio=:equal, title="hu"))
heatmap(mesh.x, mesh.y, sol.u[:,:,3], clim=clim, aspect_ratio=:equal, title="hv")

# integrator = Integrator(equation, params, Euler(), Rusanov2D(), 1, DefaultLogConfig);