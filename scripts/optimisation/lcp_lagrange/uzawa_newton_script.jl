#include("../../../src/numdiff/projection_algorithm/optimisation_process.jl")
using NumericalDiffusion
using Plots
using Test
using UnPack
using LinearAlgebra

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle
bound_mode = SingleBound()

@time sol = compute_entropic_G(params, equation; bound_mode=bound_mode, maxiter_newton=1000, use_newton=true)

@time sol = compute_entropic_G(params, equation; bound_mode=bound_mode, maxiter_uzawa=10000, use_newton=false)

# # @unpack p_uz, p_newt, w_uz, w_newt, G_uz, G_newt, D_uz, D_newt = sol
# @unpack A, b, M, q, G_newt, p_newt = sol

# # plot(p_uz, label="p uzawa only")
# # display(plot!(p_newt, label="uzawa+newton"))

# # plot(w_uz, label="w uzawa only")
# # display(plot!(w_newt, label="uzawa+newton"))

# # plot(params.mesh.x, G_uz, label="G uzawa only")
# # display(plot!(params.mesh.x, G_newt, label="uzawa+newton"))

# # plot(params.mesh.x, D_uz, label="D uzawa only")
# # display(plot!(params.mesh.x, D_newt, label="uzawa+newton"))

# # Computing residual of new solution
# @show constraint_residual = norm(max.(0.0, A*G_newt .- b))

# # Perturbating the Lagrange multiplier by staying in Ker(Aᵀ):
# c = (zero(p_newt) .+ 1) .* (p_newt.!=0) # A basis for Ker(Aᵀ)
# eps = 0.01
# ptilde = p_newt .+ eps*c

# function check_lcp(M::AbstractMatrix, q::AbstractVector, p::AbstractVector; tol=1e-10)
#     @test minimum(p) > -tol
#     @test minimum(M*p + q) > -tol
#     @test abs(p'*(M*p+q)) < tol
# end

# function show_lcp(M::AbstractMatrix, q::AbstractVector, p::AbstractVector)
#     @show minimum(p)
#     @show minimum(M*p + q)
#     @show abs(p'*(M*p+q))
# end

# show_lcp(M, q, p_newt)
# show_lcp(M, q, ptilde)

# # check_lcp(M, q, p_newt)
# # check_lcp(M, q, ptilde)

# function lcp_perturbation(M::AbstractMatrix, q::AbstractVector, p::AbstractVector, c::AbstractVector; N=10, epsmax=1.0)

#     eps = LinRange(0.0, epsmax, N)

#     # Residuals
#     pres = zeros(N)
#     wres = zeros(N)
#     complem_res = zeros(N)

#     for k in 1:N
#         #@show p
        
#         ptilde = p .+ eps[k]*c

#         #@show ptilde

#         pres[k] = minimum(min.(0.0,ptilde))
#         wres[k] = minimum(min.(0.0,M*ptilde + q))
#         complem_res[k] = (ptilde'*(M*ptilde+q))
#     end

#     eps, pres, wres, complem_res

# end

# epsmax = 0.01
# N = 10

# eps, pres, wres, complem_res = lcp_perturbation(M, q, p_newt, c; N=N, epsmax=epsmax)

# display(scatter(eps, pres, label="p residual"))
# display(scatter(eps, wres, label="w residual"))
# display(scatter(eps, complem_res, label="complementarity residual"))