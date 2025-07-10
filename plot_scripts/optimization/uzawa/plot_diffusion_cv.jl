using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
using GLM, DataFrames
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function diffusion_convergence(Nxvec::AbstractVector, alpha::Real, equation::Equation, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Real)
    N = length(Nxvec)
    diffgap = zeros(N)
    for k in 1:N
        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        # Finite volumes resolution
        sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false, false));

        # Multidimensional bounds for ΔG
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
        @unpack uinit, u, l, L, etacont_init, etacont = estimate

        # Rusanov entropic flux GRus
        GRus = G_from_theory(Rusanov(), equation, params, uinit)

        # Uzawa algorithm
        @show estimate.params.mesh.Nx
        Gc, A, b, W = init_optim_components(estimate, AlphaWeights(alpha))
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=1000000, eps=1e-8, eps_cons=1e-8);
        Gopt = optsol.gamma_opt

        #-# Numerical diffusions associated to fluxes

        DRus = zero(L)
        diffusion!(Posteriori(), GRus, etacont_init, etacont, estimate.dt, mesh, DRus)

        Dopt = zero(L)
        diffusion!(Posteriori(), Gopt, etacont_init, etacont, estimate.dt, mesh, Dopt)

        #-#

        diffgap[k] = norm(DRus .- Dopt)
    end
    
    diffgap
end

# Domain definition
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5

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

Nxvec = [30, 100, 300]
#Nxvec=[10]
alphavec = [0, 1, 2]

plt=plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

for alpha in alphavec
    diffgap = diffusion_convergence(Nxvec, alpha, equation, xmin, xmax, t0, tf, CFL_factor)
    df = DataFrame(x = log10.(Nxvec), y = log10.(diffgap))
    model = lm(@formula(y ~ x), df)
    println(model)
    coeffs = round.(coef(model); digits=2)
    scatter!(log10.(Nxvec), log10.(diffgap), label="alpha="*string(alpha)*", slope="*string(coeffs[2]), markersize=10)
    xlabel!("log(Nx)") 
end

display(title!("log(||D - Dopt||₂)"))