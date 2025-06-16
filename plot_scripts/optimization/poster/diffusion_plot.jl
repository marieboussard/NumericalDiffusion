using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
using LaTeXStrings
using GLM, DataFrames
include("../../../src/numdiff/include_file.jl")
include("../../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
alpha = 1

equation = BurgersArticle
function uexact_burgers_article(x::Real, t::Real)
    if t >= 2 / 3
        @warn "Warning: This solution is not valid for t ≥ 2/3"
    end
    if x <= -2 * t
        return -1 * (x + 2) / (1 - t)
    elseif x <= 3 * t
        return x / t
    else
        return -3 * (x - 2) / (2 - 3 * t)
    end
end

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Exact solution 
uexact = [uexact_burgers_article(xi, tf) for xi in mesh.x]

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false, false));

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L, etacont_init, etacont = estimate

# Rusanov entropic flux GRus
GRus = G_from_theory(Rusanov(), equation, params, uinit)

# Uzawa algorithm

# 1 # With two bounds
#Gc, A, b, W = init_optim_components(estimate, AbsWeights(alpha))
Gc, A, b, W = init_optim_components_upperbound_only(estimate, AbsWeights(alpha))

optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=1000000, eps=1e-8, eps_cons=1e-6);
Gopt = optsol.gamma_opt

#-# Numerical diffusions associated to fluxes

Dc = zero(L)
diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

DRus = zero(L)
diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

Dopt = zero(L)
diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

# Now, plotting Delta G and the bounds

DGRus = zero(l)
for j in 1:Nx
    DGRus[j] = (GRus[j] - GRus[mod1(j-1,Nx)])*sol.dt/mesh.dx
end

DGopt = zero(l)
for j in 1:Nx
    DGopt[j] = (Gopt[j] - Gopt[mod1(j-1,Nx)])*sol.dt/mesh.dx
end

DGc = zero(l)
for j in 1:Nx
    DGc[j] = (Gc[j] - Gc[mod1(j-1,Nx)])*sol.dt/mesh.dx
end

#-#

# Plotting

pltA=[]

plt1 = plot(size=(900, 600), margin=1Plots.cm, legend=:topleft,
legendfontsize=30,
titlefontsize=42,
guidefontsize=42,
tickfontsize=36)
plot!(mesh.x, sol.uinit, label="t=0", lw=4, lc=:black, ls=:dash)
plot!(mesh.x, sol.u, label="t = "*string(round(tf, sigdigits=2)), lw=4, lc=:red)
plot!(mesh.x, uexact, label="exact", lw=4, lc=:black)
title!("Rusanov + Euler, Nx="*string(Nx))
ylabel!("u")
push!(pltA, plt1)

plt3 = plot(size=(900, 600), margin=1Plots.cm, legend=:bottomleft,
legendfontsize=30,
titlefontsize=42,
guidefontsize=42,
tickfontsize=36)

plot!(mesh.x, l .- l, label="l - l", lw=4, lc=:black, linestyle=:dot)
plot!(mesh.x, L .- l, label="L - l", lw=4, lc=:black, linestyle=:dash)
plot!(mesh.x, DGRus .- l, label="DGRus - l", lw=4, lc=:blue)
plot!(mesh.x, DGopt .- l, label="DGopt - l", lw=4, lc=:red)
plot!(mesh.x, DGc .- l, label="DGc - l", lw=4, lc=:green)
xlabel!("x")
title!(L"\frac{Δt}{Δx}(G_{j+1/2} - G_{j-1/2})")

push!(pltA, plt3)

plt2 = plot(size=(900, 600), margin=1Plots.cm, legend=:bottomleft,
legendfontsize=30,
titlefontsize=42,
guidefontsize=42,
tickfontsize=36)
plot!(mesh.x, DRus, label="DRus", lw=4, lc=:blue)
plot!(mesh.x, Dopt, label="Dopt", lw=4, lc=:red)
plot!(mesh.x, Dc, label="Dc", lw=4, lc=:green)
#xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Numerical Diffusion")

push!(pltA, plt2)

plot(pltA..., layout=(3, 1), size=(1600, 3000))



