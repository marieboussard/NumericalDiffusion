include("../src/include_file.jl")

# # 1 # Diffusion a priori for Burgers equation

# xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
# CFL_factor = 0.5
# omega = createInterval(xmin, xmax, Nx, t0, Tf)
# eq = burgers()
# method = Rusanov(CFL_factor)

# u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)

# solEnt = optimize_for_entropy(u0, omega, eq, method)
# D_priori = diffusion_a_priori(u0, omega, eq, method)

# plot(omega.x, D_priori, label="D priori")
# plot!(omega.x, solEnt.Dopt, label="Dopt")
# xlabel!("x")
# ylabel!("Numerical Diffusion")

# 2 # Diffusion a priori for Saint Venant equation

Nx, t0, Tf = 100, 0, 0.4
CFL_factor = 0.5
domain = createUnitInterval(Nx, t0, Tf)
topography = bump_zb(height=0.5, width = 0.4)
#topography = flat_zb(height=0)
eq = SaintVenant(topography, 1e-10)
#eq = SaintVenant(NullSource(), 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
addSource!(eq.source, domain)

v0 = v0_lake_at_rest_perturbated(domain.x, eq.source)
fv_sol = fv_solve(domain, v0, eq, method)
plot_fv_sol(fv_sol, eq; nb_plots=5)

solEnt = optimize_for_entropy(v0, domain, eq, method)
D_low, D_up = diffusion_a_priori(v0, domain, eq, method)

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

plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

plot!(domain.x, D_low, label="D priori")
plot!(domain.x, solEnt.Dopt, label="Dopt")
#plot!(domain.x, D_up, label="D priori up")
xlabel!("x")
ylabel!("Numerical Diffusion")

push!(pltA, plt2)

plot(pltA..., layout=(2,1), size=(1200, 900))
