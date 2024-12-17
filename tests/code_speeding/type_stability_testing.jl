include("../src/include_file.jl")

# xmin, xmax, Nx, t0, Tf = -2, 2, 10, 0, 0.1
# CFL_factor = 0.5
# domain = createInterval(Nx, xmin, xmax, t0, Tf)
# testcase = ArticleTestcase()
# u0 = (res=zeros(domain.Nx, 1); for i in 1:Nx res[i,:]=[u0_fun(testcase, domain.x[i])] end; res)
# equation = burgers()
# method = Rusanov(CFL_factor)
# dx = domain.dx

# dt = method.CFL_factor * dx / CFL_cond(equation, u0)

# # u = zeros((6,1))
# # u[:,1] = [1,2,5,2,1,0]
# j=3
# sL, sR = 1, 1
# KFun = meanK(1,1)

# ut = compute_u_tilde(KFun, u0, j, sL, sR)

# uh = compute_u_hat(NullSource(), ut, dx, dt, j, domain, equation, method)

# m, M = initBounds(KFun, equation, u0, j, sL,sR)
# m, M = updateBounds!(KFun, NormalBounds(), equation, m, M, ut, uh, j, sL, sR, Nx, dx, dt)
# m_vec, M_vec = compute_G_bounds(u0, Nx, dx, dt, equation, domain, method, KFun)

# up = scheme_step(NullSource(), u0, dt, domain, equation, method)

# @code_warntype scheme_step(NullSource(), u0, dt, domain, equation, method)

# gamma_init = initial_guess(MeanInitGuess(), m_vec, M_vec)
# J_init = J(SquareMinFun(), gamma_init, u0, up, Nx, dx, dt, m_vec, M_vec, equation, domain)

# @code_warntype J(SquareMinFun(), gamma_init, u0, up, Nx, dx, dt, m_vec, M_vec, equation, domain)

ns, method, equation, v = NullSource(), MixedMethod{Float64, Float64}(0.5, Centered{Float64}(0.5), Rusanov{Float64}(0.5), [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0;;]), NewEq(NullSource()), [1.0; 1.0; 1.0; 1.0; 1.0; -0.75; -0.75; -0.75; -0.75; -0.75;;]
mat = giveNumFlux(ns, method, equation, v) # @enter
@code_warntype giveNumFlux(ns, method, equation, v)