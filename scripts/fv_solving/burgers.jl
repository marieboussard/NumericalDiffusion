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

sol = solve(equation, params, Euler(), Rusanov());#; log_config=LogConfig(true, true, true));

using Plots
plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
display(plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.t)))

@btime solve(equation, params, Euler(), Rusanov());


# integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig);

# #u, fu = FiniteVolumes.view_stencil!(integrator, 3)


# using FiniteVolumes: view_stencil!


# function CFL_localv2!(::Scalar, integrator::Integrator, j::Int)
#     @unpack cache, space_cache = integrator
#     @unpack cfl_cache = cache
#     @unpack absDfcont = cfl_cache
#     @unpack Nx = integrator.params.mesh

#     space_cache.cfl_loc = absDfcont[j]
#     space_cache.cfl_loc = max(space_cache.cfl_loc, absDfcont[mod1(j+1,Nx)])
# end

# @show @allocated CFL_localv2!(integrator.equation.eqtype, integrator,2)
