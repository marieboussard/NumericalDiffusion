# using Plots

# include("../src/equations/eq_saint_venant.jl")
# include("../src/domain.jl")

include("../src/include_file.jl")

eq = SaintVenant(Bump_zb())
domain = createInterval(0, 1, 100, 0, 0.4)
v0 = v0_lake_at_rest_perturbated(domain.x, eq.source)

#plot(domain.x, [v[1] for v in v0] .+ zb(eq.source, domain.x))
plot(domain.x, v0[:,1] .+ zb(eq.source, domain.x))
plot!(domain.x, zb(eq.source, domain.x))