uC_minus(a1, uR, dx) = 0.25 * (2 * a1 * dx .+ uR .- sqrt.(9 * uR .^ 2 .+ 12 * a1 * dx .* uR .+ 12 * a1^2 * dx^2))

uR_vec = collect(LinRange(0.1, 10.0, 100))
a1_vec = collect(LinRange(1, 10, 5))
dx = 0.1

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
for a1 in a1_vec
    plot!(uR_vec, uC_minus(a1, uR_vec, dx), label="a1 = " * string(a1), lw=2)
end
display(title!("Possible domain for uC"))