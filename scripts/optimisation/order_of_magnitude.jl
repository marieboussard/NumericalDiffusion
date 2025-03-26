using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function magnitude_order(Nx_vec::Vector{Int}, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Float64)
    N = length(Nx_vec)
    Gorder = zeros(Float64, N)
    dLorder = zeros(Float64, N)
    for k in 1:N
        Nx = Nx_vec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)
        equation = BurgersArticle
        # Finite volumes resolution
        sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false), maxiter=1000000);
        println(sol.status)
        # Multidimensional bounds for ΔG
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
        @unpack uinit, u, l, L = estimate
        dLorder[k] = maximum(abs.(L.-l))
        # Consistent (but non entropic) flux 
        Gc = zeros(eltype(uinit), Nx)
        Gflux!(CenteredG(), Gc, estimate)
        Gorder[k] = maximum(abs.(Gc))
        
        plot(mesh.x, 0.1*Nx^2*(L.-l), label = "L-l")
        display(plot!(mesh.x, Gc, label = "Gc", title="Nx = "*string(Nx)))
    end
    Gorder, dLorder
end

# Domain definition
xmin, xmax = -2, 2
t0, tf = 0.0, 0.2
CFL_factor = 0.5

Nx_vec = [10, 30, 100, 300, 1000, 3000, 10000, 30000]#, 100000, 300000, 1000000]
Gorder, dLorder = magnitude_order(Nx_vec, xmin, xmax, t0, tf, CFL_factor)

# x = log10.(Nx_vec)
# y = log10.(Gorder./dLorder)

# using Plots
# display(scatter(x, y))

# using GLM, DataFrames

# df = DataFrame(x = x, y = y)  # Création d'un DataFrame

# # Ajustement du modèle de régression linéaire
# model = lm(@formula(y ~ x), df)

# # Affichage des résultats
# println(model)
