using NumericalDiffusion
using LinearAlgebra
using Optim

# Test: An easy constructive algorithm in the case of total diffusion equal to zero

"""
    pA(x::AbstractVector)

Compute the orthogonal projection of the vector x onto the vector subspace Im(A)
"""
pA(x::AbstractVector) = x .- sum(x)/length(x)

"""
    check_optimality_system(gamma::AbstractVector, p::AbstractVector, A::AbstractMatrix, b::AbstractVector, f::AbstractVector; tol=1e-12)

Check that the couple (`gamma`, `p`) is solution of the optimality system parametrized by (`A`, `b`, `f`).
"""
function check_optimality_system(gamma::AbstractVector, p::AbstractVector, A::AbstractMatrix, b::AbstractVector, f::AbstractVector; tol=1e-12)

    if maximum(abs.(A'*p .+ gamma .- f)) > tol
        println("Stationarity of the lagrangian not satisfied!")
        return false
    elseif minimum(b .- A*gamma) < -tol
        
        println("Constraint Aγ-b ≤ 0 is not satisfied!")
        return false
    elseif minimum(p) < -tol
        println("Positivity of p is not satisfied!")
        return false
    elseif abs(dot(p, b.- A*gamma)) > tol
        println("Complementarity not satisfied!")
        @show b.- A*gamma
        @show p
        return false
    else
        println("Optimality system satisfied!")
        return true
    end
end

"""
    solve_optimality_system_b_in_image(b::AbstractVector, f::AbstractVector, A::AbstractMatrix; tol=1e-12)

Compute (gamma, p) solution of the optimality system
`A`^T*p + γ - `f` = 0
`A`γ - `b` ≤ 0
p ≥ 0
p^T(`A`γ-`b`) = 0

in the case where `b` is in Im(`A`). p is chosen as the solution minimal in norm.
"""
function solve_optimality_system_b_in_image(b::AbstractVector, f::AbstractVector, A::AbstractMatrix; tol=1e-12)

    m = length(b)

    # 1 # Find γ0 in Im(A) such that Aγ0 = b
    gamma_0 = pinv(A)*b
    gamma_0 = pA(gamma_0)

    # 2 # Find λ such that f - γ ∈ Ker(A), with γ = γ0 + λ*1
    lambda = dot(f.-gamma_0, ones(m))/m
    gamma = gamma_0 .+ lambda

    # 3 # Find p0 in Im(A) such that A^Tp0 = f - γ
    p0 = pinv(A')*(f-gamma)

    # 4 # Choose p = p0 + λ*1 such that p ≥ 0
    p = p0 .+ max(0.0, -minimum(p0))

    # 5 # Check that the optimality system is satisfied
    if check_optimality_system(gamma, p, A, b, f)
        println("SUCCESS")
    else
        println("FAILURE")
    end
    gamma, p
end


function solve_optimality_system_b_not_in_image(b::AbstractVector, f::AbstractVector, A::AbstractMatrix;tol=1e-12)

    m = length(b)

    # 1A # Constructs x such that x ∈ Im(A) and x ≤ b 

    # x = pA(b)

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
    # #optsol = optimize(J, ones(m-1))
    # u0 = pinv(Z)*pA(b)
    # optsol = optimize(J, u0)
    # @show optsol
    # u = Optim.minimizer(optsol)
    # @show Optim.minimum(optsol)
    # x = Z*u

    #x = [-b[2]; b[2]]

    x = zeros(m)

    @show x

    # 1B # Find γ0 in Im(A) such that Aγ0 = x
    gamma_0 = pinv(A)*x
    gamma_0 = pA(gamma_0)

    # 2 # Find λ such that f - γ ∈ Ker(A), with γ = γ0 + λ*1
    lambda = dot(f.-gamma_0, ones(m))/m
    gamma = gamma_0 .+ lambda

    # 3 # Find p0 in Im(A) such that A^Tp0 = f - γ
    p0 = pinv(A')*(f-gamma)

    # 4 # Choose p = p0 + λ*1 such that p ≥ 0
    p = p0 .+ max(0.0, -minimum(p0))

    # 5 # Check that the optimality system is satisfied
    if check_optimality_system(gamma, p, A, b, f)
        println("SUCCESS")
    else
        println("FAILURE")
    end
    gamma, p

end

# # Test data 

# #b = [1; -1]
# m = 5
# f = (rand(m) .- 1.0) .*2.0

# # Matrix A for the upper bound only
# A = zeros(m,m)
# for j in 1:m
#     A[j,j] = 1
#     for i in 1:m
#         if i == j+1
#             A[i,j] = -1
#         end
#     end    
# end
# A[1,m] = -1;

# # # Case 1: b is in Im(A)
# # b = rand(m)
# # b = pA(b)
# # gamma, p = solve_optimality_system_b_in_image(b, f, A)

# # Case 2: b is not in Im(A)
# @show b = rand(m) .+ 1.0
# gamma, p = solve_optimality_system_b_not_in_image(b, f, A);