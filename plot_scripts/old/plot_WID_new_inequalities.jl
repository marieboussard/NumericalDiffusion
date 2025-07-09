include("../src/include_file.jl");

Nx = 51
xmin, xmax = 0.0, 1.0
CFL_factor = 0.2
equation = burgers()
#cheme = FVScheme(Euler(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), Roe(CFL_factor))
scheme = FVScheme(RK2(), Rusanov(CFL_factor))
#scheme = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
#scheme = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
boxBounds = [-3 3;]
#sourceBounds=[-5.0, 5.0]
sL, sR = get_sL(scheme), get_sR(scheme)
modifiedDataType = meanK(sL, sR)
boundsType = SimpleBounds()
#boundsType = NormalBounds()

println("====== LOOKING FOR THE WORST INITIAL DATA ======")
println("Using scheme " * get_name(scheme))
@time sol = iterate_WID(xmin, xmax, Nx, equation, scheme; nb_it=1, boxBounds=boxBounds, boundsType=boundsType, modifiedDataType=modifiedDataType)
modifiedDataType = meanK(2, 2)
println("The worst value for epsilon was found to be " * string(sol.worstLowDiffVec[begin]))
#u, domain = correct_extend_initial_data(sol)
u, domain = extendInitialDataToK(Nx, sol)
dx = domain.dx
dt = scheme.spaceScheme.CFL_factor * dx / CFL_cond(equation, u)
domain.Tf = dt

println()

fv_sol = fv_solve(domain, u, equation, scheme)
#solEnt = optimize_for_entropy(u, domain, equation, method)
# method = scheme
solEnt = optimize_for_entropy(u, domain, equation, scheme; modifiedDataType=modifiedDataType, boundsType=boundsType)#scheme=LBFGS(); autodiff=:forward
println("====== APPLYING DIFFUSION QUANTIFICATION TO THIS WORST INITIAL DATA ======")
@show solEnt.optimResult

nb_plots = 2
pltA = []

# Plotting solution

sol = fv_sol
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:right,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
p = div(sol.Nt, nb_plots)
for k in 0:nb_plots-2
    plot!(sol.domain.x, sol.u_approx[k*p+1][:, 1], label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=3)
end
plot!(sol.domain.x, sol.u_approx[end][:, 1], label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=3)
xlabel!("x")
title!(get_name(sol.scheme) * ", Nx = " * string(sol.domain.Nx))
ylabel!("u")

push!(pltA, plt1)

# Numerical Entropy Flux

plt2 = plot(size=(900, 600), margin=2Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
# plot!(solEnt.domain.interfaces, solEnt.Gopt, label="Gopt", lw=2)
# plot!(solEnt.domain.interfaces, solEnt.m_vec, label="m", lw=2)

# plot!(solEnt.domain.interfaces, solEnt.M_vec, label="M", lw=2)


plot!(solEnt.domain.interfaces, solEnt.m_vec, label="m", lw=3)
plot!(solEnt.domain.interfaces, solEnt.Gopt, label="Gopt", lw=3)
plot!(solEnt.domain.interfaces, solEnt.M_vec, label="M", lw=3)
xlabel!("x")
ylabel!("Numerical Entropy Flux")
#ylims!((-1,1))

push!(pltA, plt2)

# Numerical Diffusion

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solEnt.domain.x, solEnt.Dopt, label="Dopt", lw=3)

xlabel!("x")
ylabel!("Numerical Diffusion")

push!(pltA, plt4)

# # Zoom

# plt5 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(solEnt.domain.x, solEnt.Dopt, label="Dopt", lw=2)
# xlabel!("x")
# ylabel!("Numerical Diffusion")
# title!("Zoom on positive diffusions")
# ylims!((-0.005, maximum(solEnt.Dopt)*1.2))

# push!(pltA, plt5)

# plt6 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# plot!(solEnt.domain.interfaces, solEnt.Copt, label="Copt", lw=2)
# xlabel!("x")
# ylabel!("Consistency term")

# push!(pltA, plt6)

plt5 = plot(size=(900, 600), margin=2Plots.cm, legend=:topright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(solEnt.domain.x, solEnt.Dopt, label="Dopt", lw=3)
plot!(solEnt.domain.interfaces, solEnt.Copt, label="Copt", lw=3)
xlabel!("x")
ylabel!(L"Contributions to $\mathcal{J}$")
title!("Zoom on positive diffusions")
ylims!((-0.005, max(maximum(solEnt.Dopt) * 1.2, maximum(solEnt.Copt) * 1.2)))

push!(pltA, plt5)

plot(pltA..., layout=(4, 1), size=(1400, 2300))
