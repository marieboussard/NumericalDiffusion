using NumericalDiffusion
using BenchmarkTools

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# # Burgers equation
# f(u) = u.^2/2
# Df(u) = u
# u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
# u0(x::AbstractVector) = u0.(x)
# equation = Equation(1, Scalar(), EquationFun(f, Df), u0)

equation = BurgersArticle
#time_scheme = Euler()
time_scheme = RK2()
#space_scheme = Rusanov()
space_scheme = MUSCL(Rusanov(), Minmod())

sol = solve(equation, params, time_scheme, space_scheme);#; log_config=LogConfig(true, true, true));


sol2 = solve(equation, params, time_scheme, MUSCL(Rusanov(), Superbee()))

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
fig = Figure()
ax = Axis(fig[1,1], xlabel="x", ylabel="u", title="MUSCL(Rusanov), Nx="*string(Nx))#title=get_name(space_scheme)*", Nx="*string(Nx))
lines!(ax, params.mesh.x, sol.uinit, label="t = "*string(round(sol.params.t0, sigdigits=2)))
lines!(ax, params.mesh.x, sol.u, label="Minmod, t = "*string(round(sol.t, sigdigits=2)))
lines!(ax, params.mesh.x, sol2.u, label="Superbee, t = "*string(round(sol2.t, sigdigits=2)))
lines!(ax, params.mesh.x, uexact_vec, label="exact")
axislegend(ax, position=:lt)
fig