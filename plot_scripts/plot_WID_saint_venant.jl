include("../src/include_file.jl");

Nx = 51
xmin, xmax = 0.0, 1.0
CFL_factor = 0.5
equation = SaintVenant(flat_zb(height=0.0), 1e-10)
method = createHydrostatic(CFL_factor, Rusanov)
boxBounds=[5.0 10;-5.0 5.0]
sourceBounds=[-5.0, 5.0]

@time sol = iterate_WID(xmin, xmax, Nx, equation, method; nb_it=1, boxBounds=boxBounds, sourceBounds=sourceBounds)
@show sol.worstLowDiffVec
#u, domain = correct_extend_initial_data(sol)
u, domain = extendInitialDataToK(Nx, sol)
dx = domain.dx
dt = method.CFL_factor * dx / CFL_cond(equation, u)
domain.Tf = dt

fv_sol = fv_solve(domain, u, equation, method)
#solEnt = optimize_for_entropy(u, domain, equation, method)
solEnt = optimize_for_entropy(u, domain, equation, method; method=LBFGS(), autodiff=:forward)


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
    plot!(sol.domain.x, sol.u_approx[k*p+1][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[k*p+1], sigdigits=2)), lw=3)
end
plot!(sol.domain.x, sol.u_approx[end][:,1] .+ domain.sourceVec, label="t = " * string(round(sol.t_vec[end], sigdigits=2)), lw=3)
plot!(sol.domain.x, domain.sourceVec, label="zb", lw=3)
xlabel!("x")
title!(get_name(sol.method)*", Nx = "*string(sol.domain.Nx))
ylabel!("Surface of the lake")

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

# # Zoom

# plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
#     legendfontsize=15,
#     titlefontsize=21,
#     guidefontsize=21,
#     tickfontsize=18)
# # plot!(solEnt.domain.interfaces, solEnt.Gopt, label="Gopt", lw=2)
# # plot!(solEnt.domain.interfaces, solEnt.m_vec, label="m", lw=2)

# # plot!(solEnt.domain.interfaces, solEnt.M_vec, label="M", lw=2)
# # xlabel!("x")
# # ylabel!("Numerical Entropy Flux")
# # xlims!((0.49, 0.51))
# #plot!(solEnt.domain.interfaces, solEnt.M_vec .- solEnt.m_vec, label="M-m", lw=2)

# # plot!(solEnt.domain.interfaces, (max.(0, solEnt.Gopt .- solEnt.M_vec)).^2, label="(G-M)^2", lw=2)
# # plot!(solEnt.domain.interfaces, (max.(0, solEnt.m_vec .- solEnt.Gopt)).^2, label="(m-G)^2", lw=2)
# plot!(solEnt.domain.interfaces, solEnt.Gopt .- solEnt.M_vec, label="G-M", lw=2)
# plot!(solEnt.domain.interfaces, solEnt.m_vec .- solEnt.Gopt, label="m-G", lw=2)
# ylims!((-1,1))
# # plot!(solEnt.domain.interfaces, (max.(0, solEnt.Gopt .- solEnt.M_vec)), label="(G-M)", lw=2)
# # plot!(solEnt.domain.interfaces, (max.(0, solEnt.m_vec .- solEnt.Gopt)), label="(m-G)", lw=2)

# push!(pltA, plt3)

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
ylims!((-0.005, maximum(solEnt.Dopt)*1.2))

push!(pltA, plt5)

plot(pltA..., layout=(4,1), size=(1400, 2300))
