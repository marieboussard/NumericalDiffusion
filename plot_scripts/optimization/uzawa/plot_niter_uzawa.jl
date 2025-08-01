using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
using ColorSchemes
using GLM, DataFrames
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

function niter_uzawa(Nxvec::AbstractVector, weights_type::AbstractNormWeights, equation::Equation, xmin::Real, xmax::Real, t0::Real, tf::Real, CFL_factor::Real, itermax::Int)
    N = length(Nxvec)
    niter = zeros(N)
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
        Gc, A, b, W = init_optim_components_upperbound_only(estimate, weights_type)
        optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=itermax, eps=1e-8, eps_cons=1e-6);

        niter[k] = optsol.niter

    end
    
    niter
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

Nxvec = [30, 50, 100]#, 200]
# Nxvec=[10, 30]
alphavec = [0, 0.1, 0.5, 1, 1.5]#, 1.5]

nitermat = zeros(length(alphavec), length(Nxvec))
nitermat_abs = zero(nitermat)
# coefsvec = zeros(length(alphavec))

Nalpha = length(alphavec)
itermax = 1000000

for i in 1:Nalpha
    alpha = alphavec[i]
    nitermat[i,:] = niter_uzawa(Nxvec, AlphaWeights(alpha), equation, xmin, xmax, t0, tf, CFL_factor, itermax)
    nitermat_abs[i,:] = niter_uzawa(Nxvec, AbsWeights(alpha), equation, xmin, xmax, t0, tf, CFL_factor, itermax)
    # niter = view(gapmat, i, :)
    # df = DataFrame(x = log10.(Nxvec), y = log10.(gap))
    # model = lm(@formula(y ~ x), df)
    # println(model)
    # coefsvec[i] = round.(coef(model); digits=2)[2]
     
end



# Plotting

plt=plot(size=(900, 600), margin=0.5Plots.cm, legend=:outerbottomright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
for i in 1:Nalpha
    alpha = alphavec[i]
    niter = view(nitermat, i, :)
    niter_abs = view(nitermat_abs, i, :)
    color = get(ColorSchemes.tab10.colors, i, :black)
    plot!(log10.(Nxvec), log10.(niter), color=color, label="alpha="*string(alpha), marker=:circle, markersize=10)
    @show niter_abs
    plot!(log10.(Nxvec), log10.(niter_abs), color=color, label="alpha="*string(alpha), marker=:cross, markersize=16, linestyle=:dash)
end
plot!(log10.(Nxvec), ones(length(Nxvec)).*log10(itermax), label="maxiter", lw = 3)
xlabel!("log(Nx)")
ylabel!("log(niter)")
display(title!("Number of iterations to reach convergence"))


# # The case of weights |uj+1 - uj|

# Nxvec = [30, 50, 100]#, 200]
# itermax = 1000000
# plt=plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomright,
# legendfontsize=15,
# titlefontsize=21,
# guidefontsize=21,
# tickfontsize=18)

# niter = niter_uzawa(Nxvec, AbsWeights(1), equation, xmin, xmax, t0, tf, CFL_factor, itermax)
# plot!(log10.(Nxvec), log10.(niter), label="wj=|uj+1-uj|", marker=:circle, markersize=10)

# niter = niter_uzawa(Nxvec, AlphaWeights(1), equation, xmin, xmax, t0, tf, CFL_factor, itermax)
# plot!(log10.(Nxvec), log10.(niter), label="α=1", marker=:circle, markersize=10)

# xlabel!("log(Nx)")
# ylabel!("log(niter)")
# display(title!("Number of iterations to reach convergence"))