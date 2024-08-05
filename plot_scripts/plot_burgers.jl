include("../src/include_file.jl")

# 1 # Solving Burgers equation for different choices of K

xmin, xmax, Nx, t0, Tf = -2, 2, 50, 0, 0.4
CFL_factor = 0.5
omega = createInterval(xmin, xmax, Nx, t0, Tf)
u0 = (res=zeros(omega.Nx, 1); for i in 1:Nx res[i,:]=[u0_burgers_article(omega.x[i])] end; res)
method = Rusanov(CFL_factor)
methodRoe = Roe(CFL_factor)

sL, sR = get_sL(method), get_sR(method)

#solRusanov = fv_solve(omega, u0, burgers(), method)
#plot_fv_sol(solBurgers, uexact_burgers_article)

# Different choices of K 

KVec = [meanK(sL, sR), maxK(), minK(), midLeftK(sL, sR), midRightK(sL, sR)]
labelKVec = ["Mean", "Max", "Min", "Middle Left", "Middle Right"]

solEntRusanovVec = []
solEntRoeVec = []
asymGapVec = []
asymGapRoeVec = []

# Rusanov
solAsym = optimize_for_entropy(u0, omega, burgers(), method, modifiedDataType=AsymmetricModifiedData(), method=LBFGS(), autodiff=:forward)
solAsym.label = "Asymmetric"
push!(solEntRusanovVec, solAsym)

# Roe
solAsymRoe = optimize_for_entropy(u0, omega, burgers(), methodRoe, modifiedDataType=AsymmetricModifiedData(), method=LBFGS(), autodiff=:forward)
solAsymRoe.label = "Asymmetric"
push!(solEntRoeVec, solAsymRoe)

for i in eachindex(KVec)
    # Rusanov
    solEnt = optimize_for_entropy(u0, omega, burgers(), method, modifiedDataType=KVec[i], method=LBFGS(), autodiff=:forward)
    solEnt.label = labelKVec[i]
    push!(solEntRusanovVec, solEnt)
    push!(asymGapVec, abs.((solEnt.Dopt .- solAsym.Dopt)./solAsym.Dopt))

    # Roe
    solEnt = optimize_for_entropy(u0, omega, burgers(), methodRoe, modifiedDataType=KVec[i], method=LBFGS(), autodiff=:forward)
    solEnt.label = labelKVec[i]
    push!(solEntRoeVec, solEnt)
    push!(asymGapRoeVec, abs.((solEnt.Dopt .- solAsymRoe.Dopt)./solAsymRoe.Dopt))
end

pltA = []

# Rusanov
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
sol = fv_solve(omega, u0, burgers(), method)
x = sol.domain.x
u_exact = [uexact_burgers_article(xi, sol.domain.Tf) for xi in x]
plot!(x, sol.u_approx[end], label="Approx", lw=2)
plot!(x, u_exact, label="Exact", lw=2)
xlabel!("x")
ylabel!("u")
title!("Rusanov, Nx = "*string(Nx))
#plot_fv_sol(solRusanov, uexact_burgers_article)
push!(pltA, plt1)

# Roe
plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
sol = fv_solve(omega, u0, burgers(), methodRoe)
x = sol.domain.x
u_exact = [uexact_burgers_article(xi, sol.domain.Tf) for xi in x]
plot!(x, sol.u_approx[end], label="Approx", lw=2)
plot!(x, u_exact, label="Exact", lw=2)
xlabel!("x")
ylabel!("u")
title!("Roe, Nx = "*string(Nx))
#plot_fv_sol(solRoe, uexact_burgers_article)
push!(pltA, plt2)

# Rusanov
plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottom,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
for sol in solEntRusanovVec
    plot!(omega.x, sol.Dopt, label=sol.label, lw=2)
end
xlabel!("x")
ylabel!("Numerical Diffusion")
push!(pltA, plt3)

# Roe
plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
for sol in solEntRoeVec
    plot!(omega.x, sol.Dopt, label=sol.label, lw=2)
end
xlabel!("x")
#ylabel!("Numerical Diffusion")
push!(pltA, plt4)

# Rusanov
plt5 = plot(size=(900, 600), margin=2Plots.cm, legend=:top,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
for i in eachindex(asymGapVec)
    plot!(omega.x, asymGapVec[i], label=labelKVec[i], lw=2)
end
xlabel!("x")
ylabel!("Relative Gap in Numerical Diffusion")
push!(pltA, plt5)

# Roe
plt6 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)
for i in eachindex(asymGapRoeVec)
    plot!(omega.x, asymGapRoeVec[i], label=labelKVec[i], lw=2)
end
xlabel!("x")
#ylabel!("Numerical Diffusion gap with asymmetric")
push!(pltA, plt6)

plot(pltA..., layout=(3,2), size=(1600, 1800))