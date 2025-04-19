using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function estimate_gap_gc_grus(Nxvec, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Real, equation::Equation, xc::Real)
    N = length(Nxvec)
    gap = zeros(N)
    weighted_gap = zeros(N)
    for k in 1:N
        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        # 1 # Calculating G Rusanov
        sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false))
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()))
        GRus = G_from_theory(estimate)

        # 2 # Calculating G centred
        Gc = zeros(Nx)
        Gflux!(CenteredG(), Gc, estimate)

        # 3 # Gap at index j
        j = Int(floor((xc - xmin)/mesh.dx)+1)
        gap[k] = abs(Gc[j] - GRus[j])
        weighted_gap[k] = (Gc[j] - GRus[j])^2 / (estimate.uinit[j+1]-estimate.uinit[j])^2

        # 4 # Showing position
        plot(mesh.x, estimate.uinit, label="u0")
        plot!([xc, xc], [-1, 4], label="xc")
        plot!([xmin+(j-1)*mesh.dx, xmin+(j-1)*mesh.dx], [-1, 4], label="j-1")
        plot!([xmin+j*mesh.dx, xmin+j*mesh.dx], [-1, 4], label="j")
        xlabel!("x")
        display(title!("Data visualisation"))
    end
    gap, weighted_gap
end

Nxvec = [10, 30, 100, 300, 1000, 3000, 10000]
# Domain definition
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5

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

# Position where we want to observe the convergence to a constant
xc = -3
gap, weighted_gap = estimate_gap_gc_grus(Nxvec, xmin, xmax, t0, tf, CFL_factor, equation, xc)


using GLM, DataFrames

df = DataFrame(x = log10.(Nxvec), y = log10.(gap))  # Création d'un DataFrame

# Ajustement du modèle de régression linéaire
model = lm(@formula(y ~ x), df)

# Affichage des résultats
println(model)
coeffs = round.(coef(model); digits=2)

scatter(log10.(Nxvec), log10.(gap), label ="|GRus[j+1/2] - Gc[j+1/2]|")
xlabel!("log(Nx)")
display(title!("slope="*string(coeffs[2])))

scatter!(log10.(Nxvec), log10.(weighted_gap), label ="|GRus[j+1/2] - Gc[j+1/2]|²/|u[j+1]-u[j]|²")