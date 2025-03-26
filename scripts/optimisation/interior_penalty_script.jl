using FiniteVolumes
using BenchmarkTools
using UnPack
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")
include("../../src/optimisation/interior_penalisation.jl")

# Domain definition
Nx = 20
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
W = zeros(eltype(u), Nx, Nx)
alpha=1
Gflux!(CenteredG(), Gc, estimate)
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

#eps = 0.1*Nx^2
eps=1.0

Ginit = Gexact

cache = InteriorPenCache(Nx, estimate.dt, mesh.dx, eps, W, Gc, l, L)
optsol = optimize(gamma -> JIntPen(gamma, cache), Ginit; g_tol=1e-20, iterations=20000, method=LBFGS(), autodiff=:forward)