using FiniteVolumes
using BenchmarkTools
using UnPack
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false));

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L = estimate

# Defining optimization components
Gc = zeros(eltype(u), Nx)
A = zeros(eltype(u), 2*Nx, Nx)
b = zeros(eltype(u), 2*Nx)
W = zeros(eltype(u), Nx, Nx)
alpha=1

Gflux!(CenteredG(), Gc, estimate)
fill_A!(A, estimate)
fill_b!(b, estimate)
fill_W!(W, estimate, alpha)

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


# Uzawa algorithm
optsol = optimize_uzawa(Gc, A, b; gamma0=Gexact, W=W, maxiter=100000, eps=1e-12, start_with_gamma=true);
# optsol = optimize_uzawa(Gexact, A, b; W=W, maxiter=100000, eps=1e-12);


# With posteriori estimation
estimate_opt = quantify_diffusion(sol, Posteriori(AsymmetricMD()));

# Plotting
using Plots
plot(mesh.x, estimate_opt.Gopt, label="posteriori")
plot!(mesh.x, Gc, label="consistent")
plot!(mesh.x, Gexact, label="exact")
display(plot!(mesh.x, optsol.gamma_opt, label="uzawa"))

@unpack gamma_opt = optsol
DG = zeros(Nx)
for i in eachindex(DG)
    DG[i] = (gamma_opt[i] - gamma_opt[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

plot(mesh.x, estimate.l, label="l")
plot!(mesh.x, estimate.L, label="L")
plot!(mesh.x, DGexact, label="exact")
display(plot!(mesh.x, DG, label="DG"))

# Dlow = zero(l)
D = zero(L)

@unpack etacont_init, etacont = estimate

diffusion!(Posteriori(), gamma_opt, etacont_init, etacont, estimate.dt, mesh, D)
scatter(mesh.x, estimate.D, label="priori multidim")
plot!(mesh.x, D, label="uzawa")
plot!(mesh.x, estimate_opt.D, label="posteriori")
display(title!("Max Diff:"*string(maximum(D))))

plot(optsol.popt[1:Nx], label="popt for l")
plot!(optsol.popt[Nx:end], label="popt for L")