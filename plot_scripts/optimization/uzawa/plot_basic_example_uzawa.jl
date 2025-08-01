using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 50
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

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

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

plt = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, sol.uinit, label="u0", lw=2)
xlabel!("x")
title!("Burgers with Nx = "*string(Nx))
display(plot!(mesh.x, sol.u, label="u", lw=2))

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L, etacont_init, etacont = estimate

# Rusanov entropic flux GRus
GRus = G_from_theory(Rusanov(), equation, params, uinit)

# FIRST CASE : MINIMIZING THE GAP WITH THE CENTRED FLUX #

# Uzawa algorithm
alpha = 1
Gc, A, b, W = init_optim_components(estimate, AlphaWeights(alpha))
@time optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=10000000, eps=1e-10, eps_cons=1e-10);
Gopt = optsol.gamma_opt


#-# Numerical diffusions associated to fluxes

Dc = zero(L)
diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

DRus = zero(L)
diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

Dopt = zero(L)
diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

#-#

# Plotting

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, DRus, label="DRus", lw=2)
scatter!(mesh.x, Dopt, label="Dopt", lw=2)
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
plot!(mesh.x, optsol.popt[1:Nx], label="popt for l", lw=2)
plot!(mesh.x, optsol.popt[Nx+1:end], label="popt for L", lw=2)
xlabel!("x")
ylabel!("Lagrange multiplier")
title!("Minimizing gap with centred flux")
# title!("Converged in "*string(optsol.niter)*" iterations")
push!(pltA, plt2)


# SECOND CASE : MINIMIZING THE GAP WITH THE CENTRED FLUX #

# Uzawa algotithm
@time optsol2 = optimize_uzawa(GRus, A, b; W=W, maxiter=100000, eps=1e-12, eps_cons=1e-12);
Gopt2 = optsol2.gamma_opt

#-# Numerical diffusions associated to fluxes

# Dc = zero(L)
# diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

DRus = zero(L)
diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

Dopt2 = zero(L)
diffusion!(Posteriori(), Gopt2, etacont_init, etacont, estimate.dt, mesh, Dopt2)

#-#

# Plotting

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, DRus, label="DRus", lw=2)
scatter!(mesh.x, Dopt2, label="Dopt", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Minimizing gap with entropic Rusanov flux")
push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x, optsol2.popt[1:Nx], label="popt for l", lw=2)
plot!(mesh.x, optsol2.popt[Nx+1:end], label="popt for L", lw=2)
xlabel!("x")
ylabel!("Lagrange multiplier")
title!("Minimizing gap with entropic Rusanov flux")
push!(pltA, plt4)

title_plot = plot(
    title="Quantification for Burgers with Nx = "*string(Nx)*", alpha="*string(alpha),
    titlefontsize=21,
    grid=false,
    framestyle=:none,
    xticks=false,
    yticks=false,
    # bottom_margin=-10  # pour resserrer un peu l'espace
)

# Organisation des 5 blocs (1 titre + 4 sous-plots)
layout = @layout [ a{0.05h}; [ b c; d e ] ]

# Combine le tout
plot(title_plot, pltA..., layout=layout, size=(1600, 1200))