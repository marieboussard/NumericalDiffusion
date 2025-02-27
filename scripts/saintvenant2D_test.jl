using FiniteVolumes
# using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
Ny = 100
xmin, xmax, ymin, ymax = -2, 2,-2,2
t0, tf = 0.0, 2.0
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

using FiniteVolumes: SaintVenant2Flat
equation = SaintVenant2Flat

# sol = solve(equation, params, Euler(), Rusanov2D(); maxiter=1);

# display(heatmap(mesh.x, mesh.y, sol.u[:,:,1]))
# display(heatmap(mesh.x, mesh.y, sol.u[:,:,2]))
# heatmap(mesh.x, mesh.y, sol.u[:,:,3])

integrator = Integrator(equation, params, Euler(), Rusanov2D(), 1, DefaultLogConfig);