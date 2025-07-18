using BenchmarkTools
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 50
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L = estimate

# Defining optimization components
# Gc = zeros(eltype(u), Nx)
# A = zeros(eltype(u), 2*Nx, Nx)
# b = zeros(eltype(u), 2*Nx)
# W = zeros(eltype(u), Nx, Nx)
alpha=1
bound_mode = DoubleBound()

# Gflux!(CenteredG(), Gc, estimate)
# fill_A!(A, estimate)
# fill_b!(b, estimate)
# fill_W!(W, estimate, alpha)
Gc, A, b, W = init_optim_components(bound_mode, estimate, AbsWeights(alpha))

# Exact flux
estimator = Estimator(sol, Posteriori(AsymmetricMD()), 0);
eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
eta!(estimator.entfun, estimator.u, estimator.etacont)
compute_G_bounds!(estimator)

etavec = zero(mesh.x)
Gvec = zero(mesh.x)
eta!(BurgersEnt(), estimate.uinit, etavec)
G!(BurgersEnt(), estimate.uinit, Gvec)
Gexact = zero(mesh.x)
for i in 1:Nx
    Gexact[i] = 0.5*(Gvec[i] + Gvec[mod1(i+1,Nx)]) - max(abs(estimate.uinit[i]), abs(estimate.uinit[mod1(i+1,Nx)]))*0.5*(etavec[mod1(i+1,Nx)] - etavec[i])
end
DGexact = zeros(Nx)
for i in eachindex(DGexact)
    DGexact[i] = (Gexact[i] - Gexact[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

function interior_init(gammaN::Real, params::Parameters, dt::Real, l::AbstractVector, L::AbstractVector)
    @unpack Nx, dx = params.mesh
    res = zero(l)
    alpha = sum(L)/sum(L.-l)
    res[1] = gammaN + dx/dt*(alpha*l[1]+(1-alpha)*L[1])
    for j in 2:Nx
        res[j] = res[j-1] + dx/dt*(alpha*l[j] + (1-alpha)*L[j])
    end
    res
end
gammaN = Gc[Nx]
Ginit = interior_init(gammaN, params, estimate.dt, l, L)

# Uzawa algorithm
# optsol = optimize_uzawa(Gc, A, b; gamma0=Ginit, W=W, maxiter=1000000, eps=1e-12, start_with_gamma=true);
optsol = optimize_uzawa(Gexact, A, b; W=W, maxiter=100000, eps=1e-12);
# optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=100000, eps=1e-12);


# With posteriori estimation
estimate_opt = quantify_diffusion(sol, Posteriori(AsymmetricMD()));

@unpack gamma_opt = optsol
DG = zeros(Nx)
for i in eachindex(DG)
    DG[i] = (gamma_opt[i] - gamma_opt[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

D = zero(L)
@unpack etacont_init, etacont = estimate
diffusion!(Posteriori(), gamma_opt, etacont_init, etacont, estimate.dt, mesh, D)

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

plot!(mesh.x, estimate_opt.Gopt, label="posteriori", lw=2)
plot!(mesh.x, Gc, label="consistent", lw=2)
plot!(mesh.x, Gexact, label="exact", lw=2)
plot!(mesh.x, optsol.gamma_opt, label="uzawa", lw=2)

xlabel!("x")
ylabel!("Numerical Entropy Flux")
title!("Nx="*string(Nx)*", alpha="*string(alpha))

push!(pltA, plt1)


plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottomleft,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

scatter!(mesh.x, estimate.D, label="priori multidim")
plot!(mesh.x, D, label="uzawa", lw=2)
plot!(mesh.x, estimate_opt.D, label="posteriori", lw=2)
xlabel!("x")
ylabel!("Numerical Diffusion")
title!("Max Diff:"*string(maximum(D)))
push!(pltA, plt2)


plt3 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:bottom,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

plot!(mesh.x, estimate.l, label="l", lw=2)
plot!(mesh.x, estimate.L, label="L", lw=2)
plot!(mesh.x, DGexact, label="exact", lw=2)
plot!(mesh.x, DG, label="DG", lw=2)
xlabel!("x")
ylabel!("Bounds for numerical Entropy Flux")

push!(pltA, plt3)


plt4 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)

plot!(mesh.x, optsol.popt[1:Nx], label="popt for l", lw=2)
plot!(mesh.x, optsol.popt[Nx+1:end], label="popt for L", lw=2)
xlabel!("x")
ylabel!("Lagrange multiplier")

push!(pltA, plt4)

plot(pltA..., layout=(2,2), size=(1600, 1200))