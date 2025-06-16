using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
using LaTeXStrings
using GLM, DataFrames
include("../../../src/numdiff/include_file.jl")
include("../../../src/uzawa/uzawa.jl")

function diffusion_discrepancy(Nxvec::Vector{Int}, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Real, equation::Equation)
    N = length(Nxvec)
    diff_gap = zeros(N)
    diff_tot = zeros(N)
    G_max = zeros(N)
    G_gap = zeros(N)

    for k in 1:N

        # Setting mesh size
        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        # Finite volumes resolution
        sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false, false));

        plot(mesh.x, sol.uinit)
        title!("Nx = "*string(Nx))
        display(plot!(mesh.x, sol.u))

        # Multidimensional bounds for ΔG
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
        @unpack uinit, u, l, L, etacont_init, etacont = estimate

        # Rusanov entropic flux GRus
        GRus = G_from_theory(Rusanov(), equation, params, uinit)

        # Uzawa algorithm
        # 1 # With two bounds
        Gc, A, b, W = init_optim_components(estimate, AbsWeights(1))
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=100000, eps=1e-8, eps_cons=1e-6);
        Gopt = optsol.gamma_opt

        #-# Numerical diffusions associated to fluxes

        Dc = zero(L)
        diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

        DRus = zero(L)
        diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

        Dopt = zero(L)
        diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

        # Total diffusion
        diff_tot[k] = mesh.dx*sum(abs.(DRus))

        # Discrepancy (L1 norm)
        diff_gap[k] = mesh.dx*sum(abs.(DRus .- Dopt))

        # Maximal flux 
        G_max[k] = mesh.dx*maximum(GRus)

        # Flux discrepancy (l1 norm)
        G_gap[k] = mesh.dx*sum(abs.(GRus .- Gopt))
    end

    diff_tot, diff_gap, G_max, G_gap
end


# Settings

# Domain definition
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
equation = BurgersArticle

Nxvec = [10, 30, 100, 300]

plt=plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

diff_tot, diff_gap, G_max, G_gap = diffusion_discrepancy(Nxvec, xmin, xmax, t0, tf, CFL_factor, equation)

##### NUMERICAL DIFFUSION 

# Linear regression for total diffusion
df_tot = DataFrame(x = log10.(Nxvec), y = log10.(diff_tot))
model_tot = lm(@formula(y ~ x), df_tot)
println(model_tot)
coeffs_tot = round.(coef(model_tot); digits=2)
scatter(log10.(Nxvec), log10.(diff_tot), label="total diffusion: slope="*string(coeffs_tot[2]), markersize=10, markercolor=:black)

# Linear regression for diffusion discrepancy
df_gap = DataFrame(x = log10.(Nxvec), y = log10.(diff_gap))
model_gap = lm(@formula(y ~ x), df_gap)
println(model_gap)
coeffs_gap = round.(coef(model_gap); digits=2)
scatter!(log10.(Nxvec), log10.(diff_gap), label="discrepancy: slope="*string(coeffs_gap[2]), markersize=10, markercolor=:red, marker=:cross, markerstrokewidth=4)

# Display
ylabel!("log(L1 Norm)")
title!("Rusanov + Euler")
display(xlabel!("log(Nx)"))

# ##### NUMERICAL ENTROPY FLUX 

# # Linear regression for max flux
# df_max = DataFrame(x = log10.(Nxvec), y = log10.(G_max))
# model_max = lm(@formula(y ~ x), df_max)
# println(model_max)
# coeffs_max = round.(coef(model_max); digits=2)
# scatter(log10.(Nxvec), log10.(G_max), label="Max flux: slope="*string(coeffs_max[2]), markersize=10)

# # Linear regression for flux discrepancy
# df_Ggap = DataFrame(x = log10.(Nxvec), y = log10.(G_gap))
# model_Ggap = lm(@formula(y ~ x), df_Ggap)
# println(model_Ggap)
# coeffs_Ggap = round.(coef(model_Ggap); digits=2)
# scatter!(log10.(Nxvec), log10.(G_gap), label="Discrepancy: slope="*string(coeffs_Ggap[2]), markersize=10)

# # Display
# xlabel!("log(Nx)") 
