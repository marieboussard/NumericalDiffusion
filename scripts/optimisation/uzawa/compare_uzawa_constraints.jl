using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
pyplot()
using GLM, DataFrames
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 50
# xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5
alpha = 1

# Initial condition
function u0n(n::Int, x::Real)
    if x <= -2
        return exp(-((x+3)*2)^2)/n + 1.0
    elseif x >= 2
        return exp(-((x-3)*2)^2)/n + 1.0
    elseif x <= 0
        return -2-x + 1.0
    else
        return 3-3/2*x + 1.0
    end
end
n = 2
u0(x::Real) = u0n(n, x)
u0(x::AbstractVector) = u0.(x)

# equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
xmin, xmax = -2, 2
equation = BurgersArticle

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); maxiter=10, log_config=LogConfig(true, false, true, false, false));

display(plot(mesh.x, sol.uinit))

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L, etacont_init, etacont = estimate

# Rusanov entropic flux GRus
GRus = G_from_theory(Rusanov(), equation, params, uinit)

# Uzawa algorithm
@show estimate.params.mesh.Nx

# 1 # With two bounds
Gc, A, b, W = init_optim_components(estimate, AlphaWeights(alpha))
optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=1000000, eps=1e-8, eps_cons=1e-6);
Gopt = optsol.gamma_opt

# 2 # With one bound only
Gc, A_upper, b_upper, W = init_optim_components_upperbound_only(estimate, AlphaWeights(alpha))
optsol_upper = optimize_uzawa(Gc, A_upper, b_upper; W=W, maxiter=1000000, eps=1e-8, eps_cons=1e-6);
Gopt_upper = optsol_upper.gamma_opt

#-# Numerical diffusions associated to fluxes

Dc = zero(L)
diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

DRus = zero(L)
diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

Dopt = zero(L)
diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

Dopt_upper = zero(L)
diffusion!(Posteriori(), Gopt_upper, etacont_init, etacont, estimate.dt, mesh, Dopt_upper)

#-#

# Plotting

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, DRus, label="DRus", lw=2)
plot!(mesh.x, Dopt, label="Dopt", lw=2, marker=:circle, markersize=8, markerstrokewidth = 2, markercolor = RGBA(0, 0, 0, 0), markerstrokecolor=:red, lc = :red, ls = :solid )#, markercolor = :transparent
# plot!(mesh.x, Dopt_upper, label="Dopt upper", lw=2, seriestype = :scatter, marker=(:diamond, 10), markercolor = :transparent, markerstrokecolor=:green, ls = :solid )
plot!(mesh.x, Dopt_upper, label="Dopt upper", lc=:green, lw=2, marker=:diamond, markercolor = :transparent, markerstrokecolor=:green, ls = :solid, markersize=10)
plot!(mesh.x, Dc, label="Dc", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Minimizing gap with centred flux")
push!(pltA, plt1)

plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, optsol.popt[1:Nx], lc =:red, markerstrokecolor=:red, label="popt for L", lw=2, marker=:circle, markersize=8, markerstrokewidth = 2, markercolor = :transparent, ls= :solid)
plot!(mesh.x, optsol_upper.popt[1:Nx], lc = :green, label="popt upper for L", lw=2, marker=:diamond, markersize = 10, markercolor = :transparent, markerstrokecolor=:green, ls = :solid)
# plot!(mesh.x, optsol.popt[Nx+1:end], label="popt for L", lw=2)
xlabel!("x")
ylabel!("Lagrange multiplier")
title!("Minimizing gap with centred flux")
# title!("Converged in "*string(optsol.niter)*" iterations")
push!(pltA, plt2)

title_plot = plot(
    title="Quantification for Burgers with Nx = "*string(Nx)*", alpha="*string(alpha)*", t= "*string(round(sol.t,sigdigits=3))*"s",
    titlefontsize=21,
    grid=false,
    framestyle=:none,
    xticks=false,
    yticks=false,
    # bottom_margin=-10  # pour resserrer un peu l'espace
)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:left,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
#plot!(mesh.x, DRus, label="DRus", lw=2)
plot!(mesh.x, abs.(Dopt .- DRus), markerstrokecolor=:red, lc = :red, label="|Dopt-DRus|", lw=2, marker=:circle, markersize = 8, markerstrokewidth = 2, markercolor = :transparent, ls = :solid)
plot!(mesh.x, abs.(Dopt_upper .- DRus), label="|D upper - DRus|", lw=2, marker=:diamond, markersize = 10, markercolor = :transparent, markerstrokecolor=:green, lc = :green, ls = :solid)
xlabel!("x")
ylabel!("Numerical Diffusion gap")
push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, abs.(Dopt_upper .- Dopt), markerstrokecolor=:blue, lc = :blue, label="|D upper - Dopt|", lw=2, marker=:circle, markersize = 8, markerstrokewidth = 2, markercolor = RGBA(0, 0, 0, 0), ls = :solid)
xlabel!("x")
ylabel!("Diffusion gap")
# title!("Converged in "*string(optsol.niter)*" iterations")
push!(pltA, plt4)


# Organisation des 5 blocs (1 titre + 4 sous-plots)
layout = @layout [ a{0.05h}; [ b c ; d e] ]


# Combine le tout
display(plot(title_plot, pltA..., layout=layout, size=(1600, 1200)))

plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, GRus, label="GRus", lw=2)
plot!(mesh.x, Gc, label="Gc", lw=2)
scatter!(mesh.x, Gopt, label="Gopt", marker=:circle, markersize=8)
scatter!(mesh.x, Gopt_upper, label="Gopt upper",marker=:diamond, markersize=8)
xlabel!("x")
ylabel!("Numerical entropy flux")