using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
using GLM, DataFrames
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function consistency_convergence(Nxvec::AbstractVector, alpha::Real, equation::Equation, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Real, xc::Real)
    N = length(Nxvec)
    gap = zeros(N)
    for k in 1:N
        Nx = Nxvec[k]
        mesh = OneDMesh(Nx, xmin, xmax)
        params = Parameters(mesh, t0, tf, CFL_factor)

        # Finite volumes resolution
        sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

        # Multidimensional bounds for ΔG
        estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
        @unpack uinit, u, l, L, etacont_init, etacont = estimate

        # Rusanov entropic flux GRus
        GRus = G_from_theory(Rusanov(), equation, params, uinit)

        # Uzawa algorithm
        @show estimate.params.mesh.Nx
        Gc, A, b, W = init_optim_components(estimate, AlphaWeights(alpha))
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=1000000, eps=1e-8, eps_cons=1e-6);
        Gopt = optsol.gamma_opt

        # 3 # Gap at index j
        j = Int(floor((xc - xmin)/mesh.dx)+1)
        gap[k] = abs(Gc[j] - Gopt[j])

        # 4 # Showing position
        plot(mesh.x, estimate.uinit, label="u0")
        plot!([xc, xc], [-1, 4], label="xc")
        plot!([xmin+(j-1)*mesh.dx, xmin+(j-1)*mesh.dx], [-1, 4], label="j-1")
        plot!([xmin+j*mesh.dx, xmin+j*mesh.dx], [-1, 4], label="j")
        xlabel!("x")
        display(title!("Data visualisation"))

    end
    
    gap
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
xc = -1.5

Nxvec = [30, 50, 100, 200]
# Nxvec=[10, 20]
alphavec = [0, 0.1, 1, 1.5]#, 1.5]

gapmat = zeros(length(alphavec), length(Nxvec))
coefsvec = zeros(length(alphavec))

Nalpha = length(alphavec)

for i in 1:Nalpha
    alpha = alphavec[i]
    gapmat[i,:] = consistency_convergence(Nxvec, alpha, equation, xmin, xmax, t0, tf, CFL_factor, xc)
    gap = view(gapmat, i, :)
    df = DataFrame(x = log10.(Nxvec), y = log10.(gap))
    model = lm(@formula(y ~ x), df)
    println(model)
    coefsvec[i] = round.(coef(model); digits=2)[2]
     
end



# Plotting

plt=plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

for i in 1:Nalpha
    alpha = alphavec[i]
    coeff = coefsvec[i]
    gap = view(gapmat, i, :)
    plot!(log10.(Nxvec), log10.(gap), label="alpha="*string(alpha)*", slope="*string(coeff), marker=:circle, markersize=10)
end
xlabel!("log(Nx)")
display(title!("log(|G(u(xc)) - Gopt[j]|), xc="*string(xc)))