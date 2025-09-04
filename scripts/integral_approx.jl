using Integrals
using NumericalDiffusion
using UnPack

Nx = 4
xmin, xmax = 0, 1
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
@unpack x = mesh
data = [1, 2, 1, 3]

f(t::Real) = data[min(Nx, Int(1 + floor((t * xmax + (1 - t) * xmin) / mesh.dx)))]

f(t, p) = f.(t)

# function f(t::Real)
#     @show t
#     @show Int(1 + floor((t * xmax + (1 - t) * xmin) / mesh.dx))
#     data[max(Nx, Int(1 + floor((t * xmax + (1 - t) * xmin) / mesh.dx)))]
# end

tvec = LinRange(0, 1, 10)
fvec = f(tvec, 0)

# using CairoMakie

# fig = Figure()
# ax = Axis(fig[1, 1], title="test", xlabel="t", ylabel="f")
# lines!(tvec, fvec)
# fig

function g(t::Real)
    domain = (zeros(1), t * ones(1))
    prob = IntegralProblem(f, domain)
    sol = Integrals.solve(prob, HCubatureJL(); reltol=1e-3, abstol=1e-3)
    sol.u[1]
end

g(t, p) = g.(t)

function h(x::Real)
    domain = (zeros(1), x * ones(1))
    prob = IntegralProblem(g, domain)
    sol = Integrals.solve(prob, HCubatureJL(); reltol=1e-3, abstol=1e-3)
    sol.u[1]
end

gvec = g.(tvec)
hvec = h.(tvec)

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], title="test", xlabel="t", ylabel="g")
lines!(tvec, fvec, label="f")
lines!(tvec, gvec, label="g")
lines!(tvec, hvec, label="h")
axislegend(position=:lb)
fig