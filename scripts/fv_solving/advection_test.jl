using FiniteVolumes
using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

equation = AdvectionExample

sol = solve(equation, params, Euler(), Rusanov());

plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
display(plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.params.tf)))