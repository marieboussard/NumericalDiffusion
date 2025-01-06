using Test
include("../src/tools/domain.jl")
include("../src/tools/method.jl")
include("../src/tools/equation.jl")
include("../src/tools/fv_solution.jl")
include("../src/opt_diffusion.jl")

# # Testing u tilde computation
# x = [1, 2, 3, 2, 4]
# K = 0.0
# @test compute_u_tilde(SymmetricModifiedData(), x, 2, 1, 1, K=K) == [0, 2, 3, 0, 0]
# @test compute_u_tilde(SymmetricModifiedData(), x, 5, 1, 1, K=K) == [1, 0, 0, 0, 4]
# @test compute_u_tilde(AsymmetricModifiedData(), x, 2, 1, 1) == [2, 2, 3, 3, 3]
# @test compute_u_tilde(AsymmetricModifiedData(), x, 5, 1, 1) == [1, 4, 4, 4, 4]

# # Testing u hat
# omega = createInterval(-2.0, 2.0, 10, 0.0, 0.4)
# sol = fv_solve(omega, u0_burgers_article, Burgers(), Rusanov(0.5))
# # u = sol.u_approx[end-1]
# # dt = sol.dt_vec[end]
# # dx = omega.dx
# j = 4
# sL, sR = 1, 1



# u = [0.0373636, -0.78977183, -1.10778935, -0.91451593, -0.37346463, 0.65899308,
#     1.23531824, 1.62739459, 1.47515373, 0.62020568]
# dx = 0.4
# dt = 0.12289582476537476

# ut = compute_u_tilde(meanK(sL, sR), u, j, sL, sR; K=0.0)
# uh = compute_u_hat(ut, dx, dt, j, Burgers(), Rusanov(0.5))

# # print("u:", "\n")
# # print(sol.u_approx[end-1], "\n")
# print("ut:", "\n")
# print(ut, "\n")
# print("uh:", "\n")
# print(uh, "\n")

# # Testing the cost function

CFL_number = 0.5
domain = createInterval(-2, 2, 100, 0, 0.4)
sol = optimize_for_entropy(u0_burgers_article, domain, Burgers(), Rusanov(CFL_number))
Gexact = exactG(sol.method, sol.equation, sol.u_approx[end-1])
Jexact = J(Gexact, sol.u_approx[end-1], sol.u_approx[end], sol.domain.Nx, sol.domain.dx, sol.dt_vec[end], sol.m_vec, sol.M_vec, sol.equation)

# Gexact_test = zeros(sol.domain.Nx + 1)
# Gexact_test[1:end-1] = Gexact[2:end]
# Gexact_test[end] = Gexact[1]
#Jtest = J(Gexact, sol.u_approx[end-1], sol.u_approx[end], sol.domain.Nx, sol.domain.dx, sol.dt_vec[end], sol.m_vec, sol.M_vec, sol.equation)

#@test Jexact == 0.0