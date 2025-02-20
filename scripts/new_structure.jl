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

# # Burgers equation
# f(u) = u.^2/2
# Df(u) = u
# u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
# u0(x::AbstractVector) = u0.(x)
# equation = Equation(1, Scalar(), EquationFun(f, Df), u0)

equation = BurgersArticle

@show @allocated sol = solve(equation, params, Euler(), Roe(); maxiter=2);#; log_config=LogConfig(true, true, true));

plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
display(plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.t)))


integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig)

# #u, fu = FiniteVolumes.view_stencil!(integrator, 3)


using FiniteVolumes: view_stencil!

