using BenchmarkTools
include("../../src/numdiff/include_file.jl")

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false));

estim_prio = quantify_diffusion(sol, Priori(AsymmetricMD()));
estim_priomulti = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));

using Plots
plot(mesh.x, estim_prio.D, label="priori")
display(plot!(mesh.x, estim_priomulti.D, label = "priori mulitidim"))


# Exact flux
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
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
    DGexact[i] = (Gexact[mod1(i,Nx)] - Gexact[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

plot(mesh.x, estimate.l, label="l")
plot!(mesh.x, estimate.L, label="L")
plot!(mesh.x, DGexact, label="exact")