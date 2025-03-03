using FiniteVolumes
using BenchmarkTools
using UnPack

# DOMAIN DEFINITION
Nx = 100
Ny = 100
xmin, xmax, ymin, ymax = -1, 1,-1, 1
t0, tf = 0.0, 2.0
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

# A QUADRATIC BUMP
const Az = 0.2
const Bz = 1.0
z(x, y) = max(Az - Bz*(x^2 + y^2), 0.0)
Dz(x, y) = x^2+y^2 > Az/Bz ? (zero(x), zero(x)) : (-2Bz*x, -2Bz*y)

# z(x, y) = cos(x^2 + y^2)^2 

# VISUALIZING THE TOPOGRAPHY
using Plots
x, y = mesh.x, mesh.y
c = zeros(Nx, Ny)
for j in 1:Nx
    for k in 1:Ny
        c[j,k] = z(x[j], y[k])
    end
end
display(plot(x, y, c, st=:surface, zlims=(0,1)))

# DEFINING A SAINT VENANT EQUATION
equation = saintvenant_2d_with_topo(z, Dz)

# FINITE VOLUME RESOLUTION
sol = solve(equation, params, Euler(), Rusanov2D(); maxiter=1);

clim = (minimum(sol.u .- sol.uinit), maximum(sol.u .- sol.uinit))

# display(heatmap(mesh.x, mesh.y, sol.uinit[:,:,1], clim=clim, aspect_ratio=:equal, title="h init"))
# display(heatmap(mesh.x, mesh.y, sol.uinit[:,:,2], clim=clim, aspect_ratio=:equal, title="hu init"))
# display(heatmap(mesh.x, mesh.y, sol.uinit[:,:,3], clim=clim, aspect_ratio=:equal, title="hv init"))

display(heatmap(mesh.x, mesh.y, sol.u[:,:,1].-sol.uinit[:,:,1], clim=clim, aspect_ratio=:equal, title="h - hinit"))
display(heatmap(mesh.x, mesh.y, sol.u[:,:,2], clim=clim, aspect_ratio=:equal, title="hu"))
heatmap(mesh.x, mesh.y, sol.u[:,:,3], clim=clim, aspect_ratio=:equal, title="hv")