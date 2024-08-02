include("../src/include_file.jl")

Nx = 100
CFL_factor = 0.5
eq = SaintVenant(bump_zb(height=0.5, width=0.4), 1e-10)
domain = createUnitInterval(Nx, 0.0, 0.1)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
u_init = v0_lake_at_rest(domain.x, eq.source)

# Solving Saint-Venant for this data
fv_sol = fv_solve(domain, u_init, eq, method)
display(plot_fv_sol(fv_sol, eq, nb_plots=5))

# solMeanK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, g_tol=1e-8)
# solMaxK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, modifiedDataType=maxK(), g_tol=1e-8)
# solMaxKTol = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, modifiedDataType=maxK(), g_tol=1e-16)
#solRusanov = optimize_for_entropy(u_init, domain, eq, Rusanov(CFL_factor); iterations=10000)

# Now let us try a method which is not gradient free
solMeanK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, g_tol=1e-8)
solMaxK = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-8)
solMaxKTol = optimize_for_entropy(u_init, domain, eq, method; iterations=10000, method=LBFGS(), autodiff=:forward, modifiedDataType=maxK(), g_tol=1e-10)


pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solMeanK.domain.interfaces, solMeanK.Gopt, label="Gopt", lw=2, color=:blue)
plot!(solMeanK.domain.interfaces, solMeanK.m_vec, label="m", lw=2)

plot!(solMeanK.domain.interfaces, solMeanK.M_vec, label="M", lw=2)
xlabel!("x")
title!("Optimized with Mean")
display(ylabel!("Numerical Entropy Flux"))

push!(pltA, plt1)

plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:top,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

plot!(solMaxK.domain.interfaces, solMaxK.Gopt, label="Gopt", lw=2, color=:blue)
plot!(solMaxK.domain.interfaces, solMaxK.m_vec, label="m", lw=2)
plot!(solMaxK.domain.interfaces, solMaxK.M_vec, label="M", lw=2)
xlabel!("x")
display(ylabel!("Numerical Entropy Flux"))
title!("Optimized with Max")

push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

plot!(solMeanK.domain.x, solMeanK.Dopt, label="Mean", lw=2, color=:blue)
# plot!(solMaxK.domain.x, solMaxK.Dopt, label="max K")
# plot!(solMaxKTol.domain.x, solMaxKTol.Dopt, label="max K tol")
#plot!(solRusanov.domain.x, solRusanov.Dopt, label="Rusanov")
xlabel!("x")
ylabel!("Numerical Diffusion")

push!(pltA, plt3)

# 

# plt3 = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)

# plot!(solMaxKTol.domain.interfaces, solMaxKTol.m_vec, label="m")
# plot!(solMaxKTol.domain.interfaces, solMaxKTol.Gopt, label="Gopt")
# plot!(solMaxKTol.domain.interfaces, solMaxKTol.M_vec, label="M")
# xlabel!("x")
# display(ylabel!("Numerical Entropy Flux"))

# push!(pltA, plt3)

# plt3 = plot(size=(750, 600), margin=0.5Plots.cm, legend=:bottomright,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)

# plot!(solRusanov.domain.interfaces, solRusanov.m_vec, label="m")
# plot!(solRusanov.domain.interfaces, solRusanov.Gopt, label="Gopt")
# plot!(solRusanov.domain.interfaces, solRusanov.M_vec, label="M")
# xlabel!("x")
# display(ylabel!("Numerical Entropy Flux"))

# push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:top,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

plot!(solMeanK.domain.x, solMeanK.Dopt, label="Mean (tol=1e-8)", lw=2, color=:blue)
plot!(solMaxK.domain.x, solMaxK.Dopt, label="Max (tol=1e-8)", lw=2)
plot!(solMaxKTol.domain.x, solMaxKTol.Dopt, label="Max (tol=1e-16)", lw=2)
#plot!(solRusanov.domain.x, solRusanov.Dopt, label="Rusanov")
xlabel!("x")
ylabel!("Numerical Diffusion")

push!(pltA, plt4)

plot(pltA..., layout=(2,2), size=(1600, 1200))