using NumericalDiffusion 

Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

bound_mode = SingleBound()

sol = compute_entropic_G(params, equation; bound_mode=bound_mode)

using CairoMakie 
f = Figure()
ax = Axis(f[1,1], title = get_name(sol.sol.time_scheme)*" + "*get_name(sol.sol.space_scheme), xlabel="x", ylabel="D")
lines!(ax, mesh.x, sol.D_newt, color = :tomato, label="Dopt")
scatter!(ax, mesh.x, sol.D_newt, color = :tomato)
lines!(ax, mesh.x, sol.Dc, color = :navy, label="Dc")
scatter!(ax, mesh.x, sol.Dc, color = :navy)
axislegend(position = :lb)
f