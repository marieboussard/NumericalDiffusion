using NumericalDiffusion
using BenchmarkTools
include("./beta_slopes.jl")


Nxvec = [10]
xmin, xmax = -2, 2
t0, tf = 0.0, 0.1
equation = BurgersArticle
alpha = 1
weights = AbsWeights(alpha)
bound_mode = SingleBound()

time_scheme = Euler()
space_scheme = Rusanov()

dxvec, Hvec, nitervec, maxweightvec = compute_beta_slope(Nxvec, equation; time_scheme=time_scheme, space_scheme=space_scheme)