using FiniteVolumes
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
xmin, xmax = 0, 1
t0, tf = 0.0, 0.5
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

#equation = SaintVenantAtRest
#equation = SaintVenantFlat

# SINUSOIDAL TOPOGRAPHY
const freq = 1.0
const height = 0.5
z(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
Dz(x) = pi*freq*(sin(2*pi*freq * x))*height

# # FLAT TOPOGRAPHY
# z(x) = zero(x)
# Dz(x) = zero(x)

# equation = saintvenant_with_topo(z, Dz; sourcedisc=HRDisc())
equation = saintvenant_with_topo(z, Dz)
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1);#; log_config=LogConfig(true,true,true,true));

# sol = hrsolve(params, Euler(), Rusanov(), z, Dz)

znum = z(mesh.x);

using Plots
plot(mesh.x, znum, label="topo")
display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="initial water height"))
display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

# function flux1!(y::Vector, x::Vector)
#     x[1] = y[1] + 1.0
#     x[2] = y[2] * 2.0
# end

# flux1(y::Vector) = [y[1]+1.0, y[2]*2.0]