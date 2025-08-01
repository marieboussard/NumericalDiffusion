using BenchmarkTools
include("../../src/numdiff/include_file.jl")
using FiniteVolumes:solve

# Domain definition
Nx = 6
xmin, xmax = 0, 1
t0, tf = 0.0, 0.5
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# equation = SaintVenantFlat

# # SINUSOIDAL TOPOGRAPHY
# const freq = 1.0
# const height = 0.5
# z(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
# Dz(x) = pi*freq*(sin(2*pi*freq * x))*height

# # FLAT TOPOGRAPHY
# z(x) = zero(x)
# Dz(x) = zero(x)

# QUADRATIC TOPOGRAPHY
z(x) = max.(0.0, -6*(x.-0.5).^2 .+ 0.25)
Dz(x) = -6*(x-0.5)^2+0.25 < 0.0 ? 0.0 : -12*(x-0.5)

# equation = saintvenant_with_topo(z, Dz; sourcedisc=HRDisc())
equation = saintvenant_with_topo(z, Dz)
sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false));

# sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false));

estimate = quantify_diffusion(sol, Posteriori());

znum = z(mesh.x);

using Plots
plot(mesh.x, znum, label="topo")
display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="initial water height"))
display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

plot(mesh.x, estimate.uinit, label="uinit")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))
plot(mesh.x, estimate.m, label="m")
plot!(mesh.x, estimate.M, label="M")
display(plot!(mesh.x, estimate.Gopt, label="Optimal Numerical Entropy Flux"))
plot(mesh.x, estimate.D, label="Numerical Diffusion")