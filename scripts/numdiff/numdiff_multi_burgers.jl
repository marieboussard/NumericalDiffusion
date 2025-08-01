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

estimate_priori = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()); name="priori multidim");
estimate_posteriori = quantify_diffusion(sol, Posteriori(AsymmetricMD()); name="asymmetric posteriori");

# EXACT NUMERICAL ENTROPY FLUX
estimator = Estimator(sol, Posteriori(AsymmetricMD()));
eta!(estimator.entfun, estimator.uinit, estimator.etacont_init)
eta!(estimator.entfun, estimator.u, estimator.etacont)
compute_G_bounds!(estimator)
etavec = zero(mesh.x)
Gvec = zero(mesh.x)
eta!(BurgersEnt(), estimate_posteriori.uinit, etavec)
G!(BurgersEnt(), estimate_posteriori.uinit, Gvec)
Gexact = zero(mesh.x)
for i in 1:Nx
    Gexact[i] = 0.5*(Gvec[i] + Gvec[mod1(i+1,Nx)]) - max(abs(estimate_posteriori.uinit[i]), abs(estimate_posteriori.uinit[mod1(i+1,Nx)]))*0.5*(etavec[mod1(i+1,Nx)] - etavec[i])
end
# Gback = Gexact[mod1.(Nx:1, Nx)]
Gback = [Gexact[Nx]; Gexact[1:Nx-1]]
DGexact = estimate_posteriori.dt / mesh.dx *(Gexact .- Gback)

using Plots
plot(mesh.x, estimate_priori.uinit, label="uinit")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))
plot(mesh.x, estimate_priori.l, label="l")
plot!(mesh.x, DGexact, label="DGexact")
display(plot!(mesh.x, estimate_priori.L, label="L"))
plot(mesh.x, estimate_priori.Dlow, label="Dlow")
plot!(mesh.x, estimate_posteriori.D, label=estimate_posteriori.name)
plot!(mesh.x, estimate_priori.D, label=estimate_priori.name)