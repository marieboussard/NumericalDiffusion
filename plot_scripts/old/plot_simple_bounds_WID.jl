include("../src/include_file.jl")

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
testcase = ArticleTestcase()

domain = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = initialData(domain, testcase)
scheme = FVScheme(RK2(), Rusanov(CFL_factor))
equation = burgers()
u_exact = exactData(domain, testcase)

solSimple = optimize_for_entropy(u0, domain, equation, scheme, boundsType=SimpleBounds())
sol = optimize_for_entropy(u0, domain, equation, scheme, modifiedDataType=AsymmetricModifiedData())

pltA = []

# Numerical solution
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(x, u_exact, label="Exact", lw=2)
plot!(x, solSimple.u_approx[end], label=get_name(scheme), lw=2)
xlabel!("x")
ylabel!("u")
title!(get_name(scheme) * ", Nx = " * string(Nx))
push!(pltA, plt1)

# Numerical diffusion
plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(x, solSimple.Dopt, label="Simple Bounds", lw=2)
plot!(x, sol.Dopt, label="Article Bounds", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Max Diff : " * string(maximum(solSimple.Dopt)))
push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(x, sol.Dopt, label="Article Bounds", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Max Diff : " * string(maximum(sol.Dopt)))
push!(pltA, plt3)

plot(pltA..., layout=(3, 1), size=(1000, 1200))