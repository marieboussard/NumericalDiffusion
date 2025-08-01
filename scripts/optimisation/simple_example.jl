using FiniteVolumes
using BenchmarkTools
using UnPack
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")
# include("../../src/optimisation/interior_penalty.jl")
include("interior_penalty_simple_example.jl")

N = 10
eps = 10

W = Matrix{Float64}(I, N, N)
Gc = zeros(N) .+ 2.0
b = zeros(N) .+ 1.0
A = Matrix{Float64}(I, N, N)
Ginit = zero(Gc)

cache = InteriorPenaltyCache(N, eps, W, Gc, A, b)
optsol = optimize(gamma -> J_interior_penalty(gamma, cache), Ginit; g_tol=1e-20, iterations=20000, method=LBFGS())#, autodiff=:forward)