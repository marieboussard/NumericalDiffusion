using NumericalDiffusion
using UnPack

u0(x::Real) = x<=0 ? [3.0; 0.0] : [1.0; 0.0]
# function u0(x::AbstractVector)
#     res = zeros(length(x), 2)
#     for i in eachindex(x)
#         res[i,1], res[i,2] = u0(x[i])
#     end
#     res
# end

Nx = 100
xmin, xmax = -4.0, 4.0
t0, tf = 0.0, 0.7 
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

equation = Equation(OneD(), 2, System(), SaintVenant(), u0)
time_scheme = Euler()
space_scheme = Rusanov()

sol = solve(equation, params, time_scheme, space_scheme)

@unpack x = params.mesh
hvec = [sol.u[i,1] for i in 1:Nx]
huvec = [sol.u[i,2] for i in 1:Nx]

using CairoMakie

fig = Figure(size=(1000,1000))
ax = Axis(fig[1,1], title="Water height", ylabel="h", xlabel="x")
ax2 = Axis(fig[2,1], title="Water flow", ylabel="hu", xlabel="x")
lines!(ax, x, hvec)
lines!(ax2, x, huvec)
fig