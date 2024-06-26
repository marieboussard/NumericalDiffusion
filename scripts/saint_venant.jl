using Plots

include("../src/tools/equation.jl")
include("../src/tools/method.jl")
include("../src/tools/domain.jl")
include("../src/tools/fv_solution.jl")

Nx = 100

omega = createUnitInterval(Nx, 0.0, 0.4)
#v0 = v0_lake_at_rest(omega.x, Bump_zb())
# plot(omega.x, zb(Bump_zb(), omega.x), label="zb")
# plot!(omega.x, [e[1] for e in v0] .+ zb(Bump_zb(), omega.x), label="h+zb")
# xlabel!("x")

CFL_number = 0.5
sol = fv_solve(omega, x -> v0_lake_at_rest(x, Bump_zb()), SaintVenant(Bump_zb()), Rusanov(CFL_number))

nb_plots = 5
p = div(sol.Nt, nb_plots)

# plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
#     legendfontsize=14,
#     titlefontsize=14,
#     guidefontsize=14,
#     tickfontsize=14)
plt1 = plot()

for k in 0:nb_plots
    plot!(omega.x, [sol.u_approx[k*p+1][i][1] for i in 1:Nx], label=string(round(sol.t_vec[k*p+1], sigdigits=2)))
end
xlabel!("x")
display(ylabel!("h"))

#plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright, legendfontsize=14, titlefontsize=14, guidefontsize=14, tickfontsize=14)
plt2 = plot()
for k in 0:nb_plots
    plot!(omega.x, [sol.u_approx[k*p+1][i][2] for i in 1:Nx], label=string(round(sol.t_vec[k*p+1], sigdigits=2)))
end
xlabel!("x")
display(ylabel!("hu"))

plt3 = plot()

for k in 0:nb_plots-1
    plot!(omega.x, [sol.u_approx[k*p+1][i][1] for i in 1:Nx] .+ zb(Bump_zb(), omega.x), label=string(round(sol.t_vec[k*p+1], sigdigits=2)))
end
xlabel!("x")
display(ylabel!("Surface of the lake"))