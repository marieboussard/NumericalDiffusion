include("../src/include_file.jl")

Nx = 30
CFL_factor = 0.5
eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
domain = createUnitInterval(Nx, 0.0, 0.1)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
u_init = v0_lake_at_rest(domain.x, eq.source)

println("=======NELDER MEAD======")
@time solNelderMead = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, modifiedDataType=maxK(), g_tol=1e-10)
@show solNelderMead.summary

println("======LBFGS========")
@time solLBFGS = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
@show solLBFGS.summary

println("======NEWTON========")
@time solNewton = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=Newton(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)
@show solNewton.summary

println("======A PRIORI======")
@time solPriori = diffusion_a_priori(u_init, domain, eq, method)

# plot_solution(solNelderMead)
# plot_solution(solLBFGS)
# plot_solution(solNewton)

pltA = []

plt1 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

plot!(solNelderMead.domain.x, solNelderMead.Dopt, label="Nelder Mead", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Exec time: "*string(round(solNelderMead.D, sigdigits=2)))
push!(pltA, plt1)

plt2 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solLBFGS.domain.x, solLBFGS.Dopt, label="LBFGS", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solNewton.domain.x, solNewton.Dopt, label="Newton", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solPriori.domain.x, solPriori.D_low_norm, label="Priori", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
push!(pltA, plt4)


plot(pltA..., layout=(2,2), size=(1400, 1000))