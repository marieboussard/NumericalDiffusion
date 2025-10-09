using NumericalDiffusion
using BenchmarkTools

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.2
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)


equation = BurgersArticle
time_scheme = RK2()
space_scheme = MUSCL(Rusanov(), Minmod())

# sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false));
sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false));
estimate = quantify_diffusion(sol, Posteriori(AsymmetricMD()));

sol2 = solve(equation, params, time_scheme, MUSCL(Rusanov(), Superbee()); log_config=LogConfig(true, false, true, false, false))
estimate2 = quantify_diffusion(sol2, Posteriori(AsymmetricMD()));

# Exact solution 
function uexact(x::Real, t::Real)
    if x <= -2*t
        return -(x+2)/(1-t)
    elseif x <= 3*t 
        return x/t 
    else
        return -3*(x-2)/(2-3*t)
    end
end

uexact_vec = [uexact(xi, sol.t) for xi in params.mesh.x]

using CairoMakie 
fig = Figure(size=(700,900))
ax = Axis(fig[1,1], xlabel="x", ylabel="u", title="MUSCL(Rusanov), Nx="*string(Nx))#title=get_name(space_scheme)*", Nx="*string(Nx))
#lines!(ax, params.mesh.x, sol.uinit, label="t = "*string(round(sol.params.t0, sigdigits=2)), color=:black)
lines!(ax, params.mesh.x, sol.u, label="Minmod", color=:navy)#, t = "*string(round(sol.t, sigdigits=2))
lines!(ax, params.mesh.x, sol2.u, label="Superbee", color=:tomato)
lines!(ax, params.mesh.x, uexact_vec, label="exact", color=:green)
axislegend(ax, position=:lt)

ax2 = Axis(fig[2,1], xlabel="x", ylabel="D", title="Numerical Diffusion")
lines!(ax2, params.mesh.x, estimate.D, label="Minmod", color=:navy)
lines!(ax2, params.mesh.x, estimate2.D, label="Superbee", color=:tomato)
axislegend(ax2, position=:lb)


fig