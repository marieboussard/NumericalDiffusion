using LinearAlgebra
using UnPack
using BenchmarkTools
include("../../src/uzawa/uzawa.jl")

N = 100

Gc = zeros(N) .+ 2.0
b = zeros(N) .+ 1.0
A = Matrix{Float64}(I, N, N)

sol = optimize_uzawa(Gc, A, b);