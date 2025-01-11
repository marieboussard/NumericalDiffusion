include("../src/include_file.jl")

# 1 # Solving Burgers equation for different choices of K

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
omega = createInterval(Nx, xmin, xmax, t0, Tf)
u0 = (res = zeros(omega.Nx, 1);
for i in 1:Nx
    res[i, :] = [u0_burgers_article(omega.x[i])]
end;
res)
scheme = FVScheme(Euler(), Rusanov(CFL_factor))
schemeEulerMuscl = FVScheme(Euler(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))
schemeRK2Muscl = FVScheme(RK2(), MUSCL(CFL_factor, Rusanov(CFL_factor), Minmod()))

#sL, sR = get_sL(scheme), get_sR(scheme)

#solRusanov = fv_solve(omega, u0, burgers(), scheme)
#plot_fv_sol(solBurgers, uexact_burgers_article)

#sols = []

# Rusanov
solRus = optimize_for_entropy(u0, omega, burgers(), scheme, modifiedDataType=AsymmetricModifiedData(), method=LBFGS(), autodiff=:forward)
solRus.label = "Rusanov"
#push!(sols, solRus)

# Euler + MUSCL
solEulerMuscl = optimize_for_entropy(u0, omega, burgers(), schemeEulerMuscl, modifiedDataType=AsymmetricModifiedData(), method=LBFGS(), autodiff=:forward)
solEulerMuscl.label = "Euler + MUSCL"
#push!(sols, solEulerMuscl)

# RK2 + MUSCL
solRK2Muscl = optimize_for_entropy(u0, omega, burgers(), schemeRK2Muscl, modifiedDataType=AsymmetricModifiedData(), method=LBFGS(), autodiff=:forward)
solRK2Muscl.label = "RK2 + MUSCL"
#push!(sols, solRK2Muscl)

pltA = []

# Euler MUSCL
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
sol = fv_solve(omega, u0, burgers(), scheme)
x = sol.domain.x
u_exact = [uexact_burgers_article(xi, sol.domain.Tf) for xi in x]
plot!(x, solRus.u_approx[end], label="Rusanov", lw=2)
plot!(x, solEulerMuscl.u_approx[end], label="Euler + MUSCL", lw=2)
plot!(x, u_exact, label="Exact", lw=2)
xlabel!("x")
ylabel!("u")
title!("Euler + MUSCL, Nx = " * string(Nx))
#plot_fv_sol(solRusanov, uexact_burgers_article)
push!(pltA, plt1)

# RK2 MUSCL
plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
sol = fv_solve(omega, u0, burgers(), schemeRK2Muscl)
x = sol.domain.x
u_exact = [uexact_burgers_article(xi, sol.domain.Tf) for xi in x]
plot!(x, solRus.u_approx[end], label="Rusanov", lw=2)
plot!(x, sol.u_approx[end], label="RK2 + MUSCL", lw=2)
plot!(x, u_exact, label="Exact", lw=2)
xlabel!("x")
ylabel!("u")
title!("RK2 + MUSCL, Nx = " * string(Nx))
#plot_fv_sol(solRoe, uexact_burgers_article)
push!(pltA, plt2)

# Euler MUSCL
plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottom,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(omega.x, solRus.Dopt, label=solRus.label, lw=2)
plot!(omega.x, solEulerMuscl.Dopt, label=solEulerMuscl.label, lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
ylims!(-0.0135, 0.0015)
push!(pltA, plt3)

# RK2 MUSCL
plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottom,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
plot!(omega.x, solRus.Dopt, label=solRus.label, lw=2)
plot!(omega.x, solRK2Muscl.Dopt, label=solRK2Muscl.label, lw=2)
xlabel!("x")
#ylabel!("Numerical Diffusion")
ylims!(-0.013, 0.0015)
push!(pltA, plt4)

plot(pltA..., layout=(2, 2), size=(1600, 1200))