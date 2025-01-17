include("../src/include_file.jl")

Nx = 10
CFL_factor = 0.5
a1, b1, a2, b2 = 1.0, 2.0, 2.0, 4.5
testcase = PiecewiseLinear(a1, b1, a2, b2)
xmin, xmax = spaceBounds(testcase)
@show Nx = integerNx(Nx, testcase)
equation = burgers()
scheme = FVScheme(Euler(), Roe(CFL_factor))

#1# Showing the Roe scheme is non entropic for this initial data

Tf = 0.2
domain = createInterval(Nx, xmin, xmax, 0.0, Tf)
u0 = initialData(domain, testcase)

sol = fv_solve(domain, u0, equation, scheme)
u_exact = exactData(domain, testcase)

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(domain.x, u0, label="t=0", lw=2)
plot!(domain.x, sol.u_approx[end], label=get_name(scheme) * " for t = " * string(Tf), lw=2)
plot!(domain.x, u_exact, label="Exact" * " for t = " * string(Tf), lw=2)
xlabel!("x")
ylabel!("u")
title!("Numerical resolution for Nx = " * string(Nx))
push!(pltA, plt1)

plot(pltA..., layout=(1, 1), size=(1000, 800))