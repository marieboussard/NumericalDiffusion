using FiniteVolumes
using Plots
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

# equation = Advection2Example

# SIMPLE EXAMPLE WITH KNOWN ANALYTICAL SOLUTION
# phi(v) = 1.0
# equation = advection2_vecfield(mesh, phi; sigmay=0.1)
# u0 = (x,y) -> u0_gauss2(x,y; sigmay=0.1)
# uexact = zeros(Nx, Ny)
# for j in 1:Nx
#     for k in 1:Ny
#         uexact[j,k] = exact_advection_sol(sol.t, mesh.x[j], mesh.y[k], u0)
#     end
# end

# A MORE SOPHISTICATED EXAMPLE
equation = advection2_vecfield(mesh; sigmay=0.1)

sol = solve(equation, params, Euler(), Rusanov2D())#; log_config=LogConfig(true,true,true,true));


clim = (minimum(sol.uinit), maximum(sol.uinit)) # Fixing color range
display(heatmap(mesh.x, mesh.y, sol.uinit, clim=clim, aspect_ratio=:equal, title="initial"))
# display(heatmap(mesh.x, mesh.y, uexact, clim=clim, aspect_ratio=:equal, title="exact"))
heatmap(mesh.x, mesh.y, sol.u, clim=clim, aspect_ratio=:equal, title="final")

# heatmap(mesh.x, mesh.y, sol.uinit.*a, aspect_ratio=:equal, title="fnum init")