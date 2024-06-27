using Plots

# Plot
pltA = []

zA = ["Rusanov", "Hydrostatic"]

for z in zA

    plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
        legendfontsize=14,
        titlefontsize=14,
        guidefontsize=14,
        tickfontsize=14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = " * string(z))

    xx = [1, 2, 3]
    yy = [1, 2, 3]
    yy2 = [2, 1, 1]

    scatter!(xx, yy, label="Théorique",
        markershape=:circ,
        markersize=8,
        markeralpha=nothing,
        markercolor=:red,
        markerstrokewidth=1,
        markerstrokecolor=:red
    )

    scatter!(xx, yy2, label="Numérique",
        markershape=:x,
        markersize=10,
        markeralpha=nothing,
        markercolor=:black,
        markerstrokewidth=2)

    push!(pltA, plt)

end

pltfin = plot(pltA..., layout=(1, 2), size=(1000, 500))