using NumericalDiffusion
using Plots
using LinearAlgebra
using UnPack

include("constructive_algo.jl")

Nx = 20
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle
bound_mode = SingleBound()

sol = compute_entropic_G(params, equation; bound_mode=bound_mode, maxiter_newton=1000, use_newton=true);

@unpack b, p_newt = sol

x = sol.A*sol.G_newt
x_proj = pA(b)
m = length(b)
# Z = zeros(m, m-1)
# for i in 1:m
#     for j in 1:m-1
#         if i==j
#             Z[i,j] = 1
#         elseif i==j+1 
#             Z[i,j] = -1
#         end
#     end
# end
# J(u::AbstractVector) = sum(abs.(Z*u.-b))
# u0 = pinv(Z)*pA(b)
# optsol = optimize(J, u0)
# @show optsol
# u = Optim.minimizer(optsol)
# x_l1 = Z*u


# Optimize x with the number of components equal to b 
using JuMP 
using HiGHS
using Ipopt

#model = Model(HiGHS.Optimizer)
model = Model(Ipopt.Optimizer)

@variable(model, y[i = 1:m] <= b[i])
@variable(model, z[i=1:m], Bin)
@constraint(model, sum(y[i] for i in 1:m) == 0)

# Forcer z[i] à 1 si y[i] == b[i], sinon 0 (avec tolérance ε)
ε = 1e-6
@constraint(model, [i=1:m], y[i] - b[i] >= -ε - (1 - z[i]) * 1e6)
@constraint(model, [i=1:m], y[i] - b[i] <= ε + (1 - z[i]) * 1e6)

@NLobjective(model, Max, sum(z)-sum(abs.(y.-b)))

JuMP.optimize!(model)


plot(mesh.x, b./norm(b), label="b", marker=:circle)
plot!(mesh.x, x./norm(b), label="x", marker=:circle)
plot!(mesh.x, value(y)./norm(b), label="y", marker=:circle)
# plot!(mesh.x, x_proj./norm(b), label="x proj", marker=:circle)
# plot!(mesh.x, x_l1./norm(b), label="x l1", marker=:circle)
plot!(mesh.x, p_newt./norm(p_newt), label="p", marker=:circle)