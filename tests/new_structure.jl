using FiniteVolumes
using Plots

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Burgers equation
f(u) = u.^2/2
Df(u) = u
u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
equation = Equation(1, f, Df, u0)
sol = solve(equation, params, Euler(), Rusanov())

plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.params.tf))


# integrator = FiniteVolumes.Integrator(equation, params, Euler(), Rusanov(), 100, FiniteVolumes.DefaultLogConfig)

# using UnPack
# function performstep2!(integrator::FiniteVolumes.Integrator)
#     @unpack dx = integrator.params.mesh
#     @unpack u, uprev, dt, flux = integrator
#     FiniteVolumes.numflux!(integrator)
#     @views fluxforward = flux[2:end,:]
#     @views fluxbackward = flux[1:end-1,:]
#     @. u = uprev - dt / dx * (fluxforward - fluxbackward)
# end

# using FiniteVolumes: performstep!