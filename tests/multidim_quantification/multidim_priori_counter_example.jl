include("../../src/include_file.jl")

xmin, xmax, Nx, t0 = -2, 2, 20, 0
CFL_factor = 0.5
equation = burgers()
method = Roe(CFL_factor)
#method = Rusanov(CFL_factor)
testcase = ArticleTestcase()
#testcase = SimpleShock()
domain, u0 = createOneTimestepInterval(Nx, t0, xmin, xmax, equation, testcase, CFL_factor)
modifiedDataType = minK()

D_priori = diffusion_a_priori(u0, domain, equation, method)
@show D_priori.alpha
D_priori_multidim = diffusion_a_priori_multidim(u0, domain, equation, method; modifiedDataType=modifiedDataType)
@show D_priori_multidim.alpha

m_vec, M_vec = D_priori.m_vec, D_priori.M_vec
ll_vec = domain.Tf / domain.dx * (m_vec[begin+1:end] .- M_vec[begin:end-1])
LL_vec = domain.Tf / domain.dx * (M_vec[begin+1:end] .- m_vec[begin:end-1])

l_vec, L_vec = D_priori_multidim.m_vec, D_priori_multidim.M_vec

solEnt = optimize_for_entropy(u0, domain, equation, method)

plot(domain.x, ll_vec, label="l article")
plot!(domain.x, LL_vec, label="L article")
plot!(domain.x, l_vec[begin+1:end], label="l multidim")
plot!(domain.x, L_vec[begin+1:end], label="L multidim")
display(plot!(domain.x, (solEnt.Gopt[begin+1:end] .- solEnt.Gopt[begin:end-1]) / domain.dx * domain.Tf, label="dG"))

# D_CL = D_priori.D_CL
# plot(domain.x, solEnt.Dopt, label="Dopt")
# plot!(domain.x, D_priori.D_low_norm, label="D priori low")
# plot!(domain.x, D_priori_multidim.D_low_norm, label="D priori low multidim")
# plot!(domain.x, D_priori_multidim.D_CL, label="D CL multidim")
# display(plot!(domain.x, D_CL, label="D CL"))

# plot(domain.x, solEnt.Dopt, label="Dopt")
# plot!(domain.x, D_priori.D_low, label="D low", linestyle=:dash)
# plot!(domain.x, D_priori_multidim.D_low, label="D low multidim")
# # plot!(domain.x, D_priori_multidim.D_CL, label="D CL multidim")
# plot!(domain.x, D_priori_multidim.D_low_norm, label="D norm multidim")
# plot!(domain.x, D_CL, label="D CL", linestyle=:dash)
# plot!(domain.x, D_priori.D_up, label="D up", linestyle=:dash)
# display(plot!(domain.x, D_priori_multidim.D_up, label="D up multidim"))

u_exact = [uexact_fun(testcase, xi, solEnt.domain.Tf) for xi in domain.x]
plot(domain.x, u0, label="t = " * string(t0))
plot!(domain.x, u_exact, label="u exact")
display(plot!(domain.x, solEnt.u_approx[end], label="t = " * string(domain.Tf)))

#@show uL, uN, uP = -solEnt.u_approx[1][9], -solEnt.u_approx[1][10], solEnt.u_approx[1][11]
uL, uN, uP = 1, 2, 3

@show lamb = solEnt.dt_vec[1] / domain.dx

l10 = lamb^2 / 2 * ((uP^2 - uN^2)^2 + (uN^2 - uL^2)^2) + lamb * (-uP * (uP^2 - uN^2) + uL * (uN^2 - uL^2) + uN * (uP^2 - uN^2) - uN * (uN^2 - uL^2))
@show l10