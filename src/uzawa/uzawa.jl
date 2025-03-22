using LinearAlgebra
using UnPack

include("optimizer.jl")
include("uzawa_solution.jl")
include("optimize.jl")

N = 10

Gc = zeros(N) .+ 2.0
b = zeros(N) .+ 1.0
A = Matrix{Float64}(I, N, N)

sol = optimize_uzawa(Gc, A, b)