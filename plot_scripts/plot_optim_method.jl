include("../src/include_file.jl")

Nx = 100
CFL_factor = 0.5
eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
domain = createUnitInterval(Nx, 0.0, 0.1)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
u_init = v0_lake_at_rest(domain.x, eq.source)
to = TimerOutput()

println("=======NELDER MEAD======")
@timeit to "Nelder Mead" solNelderMead = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, modifiedDataType=maxK(), g_tol=1e-10)
@show solNelderMead.summary

println("======LBFGS========")
@timeit to "LBFGS" solLBFGS = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
@show solLBFGS.summary

println("======NEWTON========")
@timeit to "Newton" solNewton = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
@show solNewton.summary

println("======A PRIORI======")
@timeit to "Priori" solPriori = diffusion_a_priori(u_init, domain, eq, method)

# plot_solution(solNelderMead)
# plot_solution(solLBFGS)
# plot_solution(solNewton)

pltA = []

plt1 = plot(size=(1200, 600), margin=1.5Plots.cm, legend=:top,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=4, ls=:dash)
plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", ls=:dot, lw=4)
plot!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=4)
plot!(solNelderMead.domain.x, solNelderMead.Dopt, label="Nelder Mead", lw=4)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))

push!(pltA, plt1)

plt2 = plot(size=(1200, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

#plot!(solNelderMead.domain.x, solNelderMead.Dopt, label="Nelder Mead", lw=2)
plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=4, ls=:dash)
plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", ls=:dot, lw=4)
plot!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=4)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))

push!(pltA, plt2)

plot(pltA..., layout=(1,2), size=(2000, 1000))

# pltA = []

# plt1 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)

# plot!(solNelderMead.domain.x, solNelderMead.Dopt, label="Nelder Mead", lw=2)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# #title!("Exec time: "*string(round(TimerOutputs.time(to["Nelder Mead"]), sigdigits=2)))
# push!(pltA, plt1)

# plt2 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=2)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# push!(pltA, plt2)

# plt3 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", lw=2)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# push!(pltA, plt3)

# plt4 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=2)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# push!(pltA, plt4)


# plot(pltA..., layout=(2,2), size=(1500, 1000))

# #### 2 #### Lake at rest perturbated

# Nx = 100
# CFL_factor = 0.5
# eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
# domain = createUnitInterval(Nx, 0.0, 0.1)
# method = createHydrostatic(CFL_factor, Rusanov)
# addSource!(eq.source, domain)
# u_init = v0_lake_at_rest_perturbated(domain.x, eq.source)
# to = TimerOutput()

# println("=======NELDER MEAD======")
# @timeit to "Nelder Mead" solNelderMead = optimize_for_entropy(u_init, domain, eq, method; iterations=10000)
# @show solNelderMead.summary

# println("======LBFGS========")
# @timeit to "LBFGS" solLBFGS = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward)
# @show solLBFGS.summary

# println("======NEWTON========")
# @timeit to "Newton" solNewton = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward)
# @show solNewton.summary

# println("======A PRIORI======")
# @timeit to "Priori" solPriori = diffusion_a_priori(u_init, domain, eq, method)

# # plot_solution(solNelderMead)
# # plot_solution(solLBFGS)
# # plot_solution(solNewton)

# pltA = []

# plt1 = plot(size=(1200, 600), margin=1.5Plots.cm, legend=:bottomleft,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)

# plot!(solNelderMead.domain.x, solNelderMead.Dopt, label="Nelder Mead", lw=2)
# plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=2, ls=:dash)
# plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", ls=:dot, lw=2)
# scatter!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=2, ms=4)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))
# #title!("Exec time: "*string(round(TimerOutputs.time(to["Nelder Mead"]), sigdigits=2)))
# #push!(pltA, plt1)

# # plt2 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
# #     legendfontsize=15,
# #     titlefontsize=21,
# #     guidefontsize=21,
# #     tickfontsize=18)
# # plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=2)
# # xlabel!("x")
# # ylabel!("Numerical Diffusion")
# # push!(pltA, plt2)

# # plt3 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
# #     legendfontsize=15,
# #     titlefontsize=21,
# #     guidefontsize=21,
# #     tickfontsize=18)
# # plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", lw=2)
# # xlabel!("x")
# # ylabel!("Numerical Diffusion")
# # push!(pltA, plt3)

# # plt4 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
# #     legendfontsize=15,
# #     titlefontsize=21,
# #     guidefontsize=21,
# #     tickfontsize=18)
# # plot!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=2)
# # xlabel!("x")
# # ylabel!("Numerical Diffusion")
# # push!(pltA, plt4)


# #plot(pltA..., layout=(2,2), size=(1500, 1000))