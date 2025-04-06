using FiniteVolumes
using BenchmarkTools
using UnPack

# Domain definition
Nx = 1000
xmin, xmax = -2, 3
t0, tf = 0.0, 1.0
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

u0(n::Int, x::Real) = exp(-x^2)/n

using Plots
plt=plot()
for n in 1:10
    plot!(mesh.x, [u0(n,xi) for xi in mesh.x],label="n="*string(n))
end
display(plt)

# function u0_test(x::Real)
#     if x <= 0
#         return zero(x)
#     elseif x >= 1
#         return zero(x)
#     else
#         return one(x)*1.0
#     end
# end

# function u0_test(x::Real)
#     if x<=0
#         return zero(x)
#     elseif x >=1
#         return zero(x)
#     else
#         return x
#     end
# end

# function u0_test(x::Real)
#     n=1
#     return exp(-n*x^2)/n
# end
# u0_test(x::AbstractVector) = u0_test.(x)
# equation = Equation(OneD(), 1, Scalar(), Burgers(), u0_test)

# sol = solve(equation, params, Euler(), Rusanov());#; log_config=LogConfig(true, true, true));

# using Plots
# plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
# display(plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.t)))