using FiniteVolumes
using BenchmarkTools
using UnPack
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")
include("../../src/optimisation/interior_penalty.jl")

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

eps = 1/(0.1*Nx^2)
# eps=1.0
gammaN = Gc[Nx]

Ginit = interior_init(gammaN, params, estimate.dt, l, L)

DG = zeros(Nx)
for i in eachindex(DG)
    DG[i] = (Ginit[i] - Ginit[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

using Plots
plot(mesh.x, Gc, label="Gc")
display(plot!(mesh.x, Ginit, label="Ginit"))

plot(mesh.x, l, label="l")
plot!(mesh.x, DG, label="DG")
display(plot!(mesh.x, L, label="L"))

cache = InteriorPenCache(Nx, estimate.dt, mesh.dx, eps, W, Gc, l, L)
@show JIntPen(Ginit, cache)
optsol = optimize(gamma -> JIntPen(gamma, cache), Ginit; g_tol=1e-20, iterations=20000, method=LBFGS(), autodiff=:forward)