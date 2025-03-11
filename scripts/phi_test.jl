using FiniteVolumes
using BenchmarkTools
using UnPack

const a = 1
const alpha = sqrt(1+a/2)

# Domain definition
Nx = 100
Ny = 100
xmin, xmax, ymin, ymax = -(1+a), 1+a, -(1+a), (1+a)
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

function P(v::Real)
    c = 1/((1-alpha)^2*(1/3+1/2*(alpha/3-1)))
    c*(v-alpha)^2*(1/3*v+1/2*(alpha/3-1))
end

function phi(v::Real)
    if v ≤ 1
        return 1
    elseif v ≥ alpha
        return 0
    else
        return P(v)
    end
end
# v = collect(LinRange(0, sqrt(1+a), 1000))
# ph = phi.(v)
# using Plots
# plot(v, ph)
# plot!([alpha, alpha], [0, 1])

equation = advection2_vecfield(mesh, phi; sigmax=0.7, sigmay=0.3)
sol = solve(equation, params, Euler(), Rusanov2D())

using Plots
clim = (minimum(sol.uinit), maximum(sol.uinit)) # Fixing color range
display(heatmap(mesh.x, mesh.y, sol.uinit, clim=clim, aspect_ratio=:equal, title="initial"))
heatmap(mesh.x, mesh.y, sol.u, clim=clim, aspect_ratio=:equal, title="final")