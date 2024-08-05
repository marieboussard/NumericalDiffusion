include("../src/include_file.jl")
Nx, t0, Tf = 100, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
topography = bump_zb(height=0.5, width = 0.4)
eq = SaintVenant(topography, 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)
v0 = v0_lake_at_rest_perturbated(domain.x, eq.source)

fv_sol = fv_solve(domain, v0, eq, method)
plot_fv_sol(fv_sol, eq; nb_plots=5)

solEnt = optimize_for_entropy(v0, domain, eq, method, method=LBFGS(), autodiff=:forward)
D_priori = diffusion_a_priori(v0, domain, eq, method)

D_low = D_priori.D_low_norm
D_up = D_priori.D_up_norm

nb_plots = 5
pltA = []

sol = fv_sol
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:right,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
p = div(sol.Nt, nb_plots)
for k in 0:nb_plots-2
    plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=2)
end
plot!(sol.domain.x, sol.u_approx[end][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=2)
plot!(sol.domain.x, domain.sourceVec, label="zb", lw=2)
xlabel!("x")
title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))
ylabel!("Surface of the lake")

push!(pltA, plt1)

plt2 = plot(size=(900, 600), margin=1.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)


plot!(domain.x, solEnt.Dopt, label="Dopt", lw=2)
scatter!(domain.x, D_low, label="D priori")#, markershape = :circs, lw=2)
#plot!(domain.x, D_up, label="D priori up")
xlabel!("x")
ylabel!("Numerical Diffusion")

push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=18,
    tickfontsize=18)


#plot!(domain.x, abs.(solEnt.Dopt, label="Dopt", lw=2)
scatter!(domain.x, abs.((solEnt.Dopt .- D_low)./solEnt.Dopt), label=L"$\frac{|Dopt-Dpriori|}{|Dopt|}$")#, markershape = :circ, lw=2)
#plot!(domain.x, D_up, label="D priori up")
xlabel!("x")
ylabel!("Discrepancy between quantifications")

push!(pltA, plt3)

plot(pltA..., layout=(3,1), size=(1200, 1400))
