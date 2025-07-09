include("../src/include_file.jl")

# Domain
xmin, xmax, Nx, t0, Tf = -2, 2, 100, 0, 0.25
CFL_factor = 0.5
domain = createInterval(Nx, xmin, xmax, t0, Tf)
equation = burgers()
scheme = FVScheme(Euler(), Roe(CFL_factor))
#scheme = FVScheme(Euler(), Rusanov(CFL_factor))

testcase = ArticleTestcase()
u0 = initialData(domain, testcase)

FVsol = fv_solve(domain, u0, equation, scheme)
u_approx, dt_vec = FVsol.u_approx, FVsol.dt_vec
dx = domain.dx
u_exact = [uexact_fun(testcase, xi, Tf) for xi in domain.x]

#modifiedDataType = midLeftK_multidim(1,1)
#modifiedDataType = minK()
#modifiedDataType = meanK_multidim(1,1)
modifiedDataType = AsymmetricModifiedData()

priori_sol = diffusion_a_priori(u0, domain, equation, scheme; modifiedDataType=modifiedDataType)
priori_sol_multidim = diffusion_a_priori_multidim(u0, domain, equation, scheme; modifiedDataType=modifiedDataType)

m_vec, M_vec = priori_sol.m_vec, priori_sol.M_vec
ll_vec = dt_vec[end] / dx * (m_vec[begin+1:end] .- M_vec[begin:end-1])
LL_vec = dt_vec[end] / dx * (M_vec[begin+1:end] .- m_vec[begin:end-1])

l_vec, L_vec = priori_sol_multidim.l_vec, priori_sol_multidim.L_vec

solEnt = optimize_for_entropy(u0, domain, equation, scheme)

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(domain.x, u0, label="t = " * string(t0), lw=2, color=:blue)
plot!(domain.x, u_exact, label="u exact", lw=2, color=:red)
plot!(domain.x, solEnt.u_approx[end], label="t = " * string(Tf), lw=2, color=:green)
xlabel!("x")
ylabel!("u")
title!("FV resolution with scheme " * get_name(scheme))
push!(pltA, plt1)

# plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(domain.x, ll_vec, label="l article", lw=2)
# plot!(domain.x, LL_vec, label="L article", lw=2)
# plot!(domain.x, l_vec[begin+1:end], label="l multidim", lw=2)
# plot!(domain.x, L_vec[begin+1:end], label="L multidim", lw=2)
# plot!(domain.x, (solEnt.Gopt[begin+1:end] .- solEnt.Gopt[begin:end-1]) / domain.dx * dt_vec[end], label="dG", lw=2)
# xlabel!("x")
# push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(domain.x, solEnt.Dopt, label="Dopt", lw=2, color=:blue)
plot!(domain.x, priori_sol_multidim.D_low, label="D low multidim", color=:purple, lw=2)
plot!(domain.x, priori_sol.D_low, label="D low", linestyle=:dash, color=:orange, lw=2)
plot!(domain.x, priori_sol_multidim.D_up, label="D up multidim", color=:pink, lw=2)
plot!(domain.x, priori_sol.D_up, label="D up", linestyle=:dash, color=:skyblue, lw=2)
plot!(domain.x, priori_sol_multidim.D_priori, label="D priori multidim", color=:red, lw=2)
plot!(domain.x, priori_sol.D_low_norm, linestyle=:dash, label="D priori", color=:green, lw=2)
# plot!(domain.x, D_priori_multidim.D_CL, label="D CL multidim")
#plot!(domain.x, D_CL, label="D CL", linestyle=:dash, lw=2)
xlabel!("x")
title!("Bounds on diffusion and diffusion estimates")
push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(domain.x, solEnt.Dopt, label="Dopt", lw=2, color=:blue)
plot!(domain.x, priori_sol_multidim.D_priori, label="D priori multidim", lw=2, color=:red)
plot!(domain.x, priori_sol.D_low_norm, linestyle=:dash, label="D priori", lw=2, color=:green)
xlabel!("x")
title!("Zoom on diffusion estimates")
push!(pltA, plt4)

plot(pltA..., layout=(3, 1), size=(1000, 1400))