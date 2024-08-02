include("../src/include_file.jl")

Nx = 30
CFL_factor = 0.5
eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
domain = createUnitInterval(Nx, 0.0, 0.1)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
u_init = v0_lake_at_rest(domain.x, eq.source)

# Solving Saint-Venant for this data
fv_sol_hydro = fv_solve(domain, u_init, eq, method)
fv_sol_rus = fv_solve(domain, u_init, eq, Rusanov(CFL_factor))

println("=======NELDER MEAD======")
@time solMeanK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000)
@show solMeanK.summary
@time solRusanov = optimize_for_entropy(u_init, domain, eq, Rusanov(CFL_factor); iterations=10000, g_tol=1e-10)

println("======LBFGS========")
@time solMeanK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, g_tol=1e-10)
@show solMeanK.summary
# solRusanov = optimize_for_entropy(u_init, domain, eq, Rusanov(CFL_factor); iterations=10000, method=LBFGS(), autodiff=:forward)
println("======NEWTON========")
@time solMeanK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward, g_tol=1e-10);
@show solMeanK.summary
#solRusanov = optimize_for_entropy(u_init, domain, eq, Rusanov(CFL_factor); iterations=10000, method=Newton(), autodiff=:forward)

# nb_plots = 5
# pltA = []

# sol = fv_sol_hydro
# plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:right,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# p = div(sol.Nt, nb_plots)
# for k in 0:nb_plots-2
#     plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
# end
# plot!(sol.domain.x, sol.u_approx[end][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
# plot!(sol.domain.x, domain.sourceVec, label="zb", lw=2)
# xlabel!("x")
# title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))
# ylabel!("Surface of the lake")

# push!(pltA, plt1)

# sol = fv_sol_rus
# plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:right,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)
# p = div(sol.Nt, nb_plots)
# for k in 0:nb_plots-2
#     plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
# end
# plot!(sol.domain.x, sol.u_approx[end][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
# plot!(sol.domain.x, domain.sourceVec, label="zb", lw=2)
# xlabel!("x")
# title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))
# #ylabel!("Surface of the lake")

# push!(pltA, plt2)

# sol=fv_sol_hydro
# plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:top,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)
# p = div(sol.Nt, nb_plots)
# for k in 0:nb_plots-2
#     plot!(sol.domain.x, sol.u_approx[k*p+1][:,2], label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
# end
# plot!(sol.domain.x, sol.u_approx[end][:,2], label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
# xlabel!("x")
# ylabel!("Water flow")
# ylims!((-0.12, 0.12))
# push!(pltA, plt3)

# sol=fv_sol_rus
# plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:top,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)
# p = div(sol.Nt, nb_plots)
# for k in 0:nb_plots-2
#     plot!(sol.domain.x, sol.u_approx[k*p+1][:,2], label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
# end
# plot!(sol.domain.x, sol.u_approx[end][:,2], label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
# xlabel!("x")
# #ylabel!("Water flow")
# ylims!((-0.12, 0.12))
# push!(pltA, plt4)

# plt5 = plot(size=(900, 600), margin=2Plots.cm, legend=:top,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)
# plot!(solMeanK.domain.x, solMeanK.Dopt, label="Dopt", lw=2)
# # plot!(solMaxK.domain.x, solMaxK.Dopt, label="max K")
# # plot!(solMaxKTol.domain.x, solMaxKTol.Dopt, label="max K tol")
# #plot!(solRusanov.domain.x, solRusanov.Dopt, label="Rusanov")
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# ylims!((-0.015, 0.009))
# push!(pltA, plt5)


# plt6 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:top,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)

# plot!(solRusanov.domain.x, solRusanov.Dopt, label="Dopt", lw=2)
# xlabel!("x")
# #ylabel!("Numerical Diffusion")
# ylims!((-0.015, 0.009))
# push!(pltA, plt6)

# plot(pltA..., layout=(3,2), size=(1600, 1800))