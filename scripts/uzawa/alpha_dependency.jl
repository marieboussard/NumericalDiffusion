using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function alpha_dependency(alpha_vec::AbstractVector, params::Parameters, equation::Equation)
    # Finite volumes resolution
    sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));
    plot(mesh.x, sol.uinit, label="u0")
    display(plot!(mesh.x, sol.u, label="u"))

    # Multidimensional bounds for ΔG
    estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
    @unpack uinit, u, l, L, etacont_init, etacont = estimate

    


    # Defining optimization components
    Gc = zeros(eltype(u), Nx)
    A = zeros(eltype(u), 2*Nx, Nx)
    b = zeros(eltype(u), 2*Nx)
    W = zeros(eltype(u), Nx, Nx)

    Gflux!(CenteredG(), Gc, estimate)
    fill_A!(A, estimate)
    fill_b!(b, estimate)

    # Diffusion of Gc
    Dc = zero(L)
    diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

    # Exact flux
    Gexact = G_from_theory(Rusanov(), equation, params, uinit)
    DGexact = zeros(Nx)
    for i in eachindex(DGexact)
        DGexact[i] = (Gexact[i] - Gexact[mod1(i-1,Nx)])*estimate.dt/mesh.dx
    end

    # Diffusion for exact flux
    Dexact = zero(L)
    diffusion!(Posteriori(), Gexact, etacont_init, etacont, estimate.dt, mesh, Dexact)

    # Pointwise evaluation of G
    Gcont = zero(uinit)
    G!(entropy(equation.funcs), uinit, Gcont)

    N = length(alpha_vec)
    mu_vec = zeros(N)
    consistency_gap = zeros(N)
    constraint_res = zeros(N)
    maxdiff = zeros(N)
    D = zero(L)

    for k in 1:N
        alpha = alpha_vec[k]
        W = zeros(eltype(u), Nx, Nx)
        fill_W!(W, estimate, alpha)
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=500000, eps=1e-12, eps_cons=1e-12);
        mu_vec[k] = optsol.mu
        window = 1:10
        consistency_gap[k] = norm(optsol.gamma_opt[window] .- Gcont[window])
        constraint_res[k] = optsol.constraint_residual
        diffusion!(Posteriori(), optsol.gamma_opt, etacont_init, etacont, estimate.dt, mesh, D)
        maxdiff[k] = maximum(D)

        plot(mesh.x[window], optsol.gamma_opt[window], label="uzawa")
        scatter!(mesh.x[window], Gexact[window], label="from theory", marker=:square)
        scatter!(mesh.x[window], Gc[window], label="centred", marker=:diamond)
        scatter!(mesh.x[window], Gcont[window], label="G(u)", marker=:cross)
        xlabel!("x")
        ylabel!("G")
        display(title!("Entropy flux in the constant area for α = "*string(alpha)))
    end
    mu_vec, consistency_gap, constraint_res, maxdiff
end

# Domain definition
Nx = 100
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersConstant

alpha_vec = [0, 0.01, 0.1, 0.5, 1, 2]
# alpha_vec = [0]
mu_vec, consistency_gap, constraint_res, maxdiff = alpha_dependency(alpha_vec, params, equation)

# Visualisation
pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
scatter!(alpha_vec, mu_vec, label="μ", marker=(:circle, 8))
xlabel!("α")
push!(pltA, plt1)

plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
scatter!(alpha_vec, log10.(consistency_gap), label="log||Gopt-G(u)||₂", marker=(:circle, 8))
xlabel!("α")
push!(pltA, plt2)

plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
scatter!(alpha_vec, log10.(constraint_res), label="log||max(0, Aγ - b)||₂", marker=(:circle, 8))
xlabel!("α")
push!(pltA, plt3)

plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
scatter!(alpha_vec, log10.(maxdiff), label="log(max(Dopt))", marker=(:circle, 8))
xlabel!("α")
push!(pltA, plt4)

# plttot = plot(pltA..., layout=(2,2), size=(1600, 1200))
# annotate!(plttot, 0.5, 1.05, text("Quantification for Burgers with Nx = "*string(Nx), :center, 16, :black))

title_plot = plot(
    title="Quantification for Burgers with Nx = "*string(Nx),
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