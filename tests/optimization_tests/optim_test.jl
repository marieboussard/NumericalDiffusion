#using Optimization
using Optim
# using OptimizationBBO


fcost(x) = sum(x .^ 2)
x_init = ones(5)
#x_init = [5.0]

@show sol = optimize(x -> fcost(x), x_init; g_tol=1e-16, method=LBFGS(), autodiff=:forward)

@show Optim.minimizer(sol)