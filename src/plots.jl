# using Plots

# include("../src/tools/fv_solution.jl")

include("include_file.jl")

Nx = 20

omega = createUnitInterval(Nx, 0.0, 0.4)
v0 = v0_lake_at_rest(omega.x, Bump_zb())
v0_perturbated = v0_lake_at_rest_perturbated(omega.x, Bump_zb(); height=0.1)

CFL_number = 0.5

# Classic Rusanov Method
sol_Rus = fv_solve(omega, v0, SaintVenant(Bump_zb()), Rusanov(CFL_number))
sol_Rus_perturbated = fv_solve(omega, v0_perturbated, SaintVenant(Bump_zb()), Rusanov(CFL_number))

# Hydrostatic Method
sol_Hydro = fv_solve(omega, v0, SaintVenant(Bump_zb()), Hydrostatic(CFL_number, Rusanov(CFL_number)))
sol_Hydro_perturbated = fv_solve(omega, v0_perturbated, SaintVenant(Bump_zb()), Hydrostatic(CFL_number, Rusanov(CFL_number)))

# Plot

nb_plots = 5


pltA = []
label_vec = ["Rusanov", "Hydrostatic", "Rusanov Perturbated", "Hydrostatic Perturbated"]
sol_vec = [sol_Rus, sol_Hydro, sol_Rus_perturbated, sol_Hydro_perturbated]

for j in eachindex(sol_vec)

    p = div(sol_vec[j].Nt, nb_plots)

    plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
        legendfontsize=10,
        titlefontsize=14,
        guidefontsize=14,
        tickfontsize=14)

    xlabel!("x")
    ylabel!("Water surface")
    title!(label_vec[j] * " " * string(Nx))

    for k in 0:nb_plots-1
        plot!(omega.x, [sol_vec[j].u_approx[k*p+1][i][1] for i in 1:Nx] .+ zb(Bump_zb(), omega.x), label="t = " * string(round(sol_vec[j].t_vec[k*p+1], sigdigits=2)))
    end
    plot!(omega.x, zb(Bump_zb(), omega.x), label="zb")

    push!(pltA, plt)

end

pltfin = plot(pltA..., layout=(2, 2), size=(1200, 1000))


# # Water flow

# pltB = []

# for j in eachindex(sol_vec)

#     p = div(sol_vec[j].Nt, nb_plots)

#     plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
#         legendfontsize=10,
#         titlefontsize=14,
#         guidefontsize=14,
#         tickfontsize=14)

#     xlabel!("x")
#     ylabel!("Water flow")
#     title!(label_vec[j] * " " * string(Nx))

#     for k in 0:nb_plots-1
#         plot!(omega.x, [sol_vec[j].u_approx[k*p+1][i][2] for i in 1:Nx], label="t = " * string(round(sol_vec[j].t_vec[k*p+1], sigdigits=2)))
#     end

#     push!(pltB, plt)

# end

# pltfinB = plot(pltB..., layout=(2, 2), size=(1200, 1000))