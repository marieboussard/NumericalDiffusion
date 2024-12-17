# using Logging, BenchmarkTools#, #ECOSProblems, JSOSolvers, NLPModels
# using BoundedLeastSquares # from https://github.com/nboyd/BoundedLeastSquares.jl

# global_logger(NullLogger())

# m, n = 10000, 10
# A = randn(m, n); x = randn(n); b = A*x;
# bl = zeros(n); bu = ones(n)
# Q = form_quadratic_from_least_squares(A, b)
# K = sqrt(Q.Q); b_K = K\(Q.b)
# x_opt = min_bound_constrained_quadratic(Q,bl,bu)

using JuMP, OSQP

model = Model(OSQP.Optimizer)
@variable(model, x)
@variable(model, y)
@objective(model, Min, x^2+y^2)
#@constraint(model, c1, x >= y+1)

print(model)
optimize!(model)