using FiniteVolumes
using BenchmarkTools
using UnPack
using Plots
include("../../src/numdiff/include_file.jl")
include("../../src/uzawa/uzawa.jl")

# Domain definition
Nx = 100
xmin, xmax = -4, 4
t0, tf = 0.0, 0.5
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersConstant
# u0(x::Real) = exp(-2*x^2)
# u0(x::AbstractVector) = u0.(x)
# # equation = Equation(OneD(), 1, Scalar(), Advection(2), u0)
# equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, false, true, false, false));

# plot(mesh.x, sol.uinit, label="u0")
# display(plot!(mesh.x, sol.u, label="u"))

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()));
@unpack uinit, u, l, L, etacont, etacont_init = estimate

# Defining optimization components
Gc = zeros(eltype(u), Nx)
A = zeros(eltype(u), 2*Nx, Nx)
b = zeros(eltype(u), 2*Nx)
W = zeros(eltype(u), Nx, Nx)
alpha=1.0

Gflux!(CenteredG(), Gc, estimate)
fill_A!(A, estimate)
fill_b!(b, estimate)
fill_W!(W, estimate, alpha)
display(plot(mesh.x, [W[i,i] for i in 1:Nx], label="w"))

# Exact flux
Gexact = G_from_theory(Rusanov(), equation, params, uinit)
DGexact = zeros(Nx)
for i in eachindex(DGexact)
    DGexact[i] = (Gexact[i] - Gexact[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

# Diffusion for exact flux
Dexact = zero(L)
diffusion!(Posteriori(), Gexact, etacont_init, etacont, estimate.dt, mesh, Dexact)

# Pointwise evaluation of G
Gcont = zero(uinit)
G!(entropy(equation.funcs), uinit, Gcont)

# Uzawa algorithm
# optsol = optimize_uzawa(Gexact, A, b; W=W, maxiter=100000, eps=1e-12);
optsol = optimize_uzawa(Gc, A, b; W=W, maxiter=500000, eps=1e-12, eps_cons=1e-12);
@show optsol.mu
@unpack gamma_opt = optsol
DG = zeros(Nx)
for i in eachindex(DG)
    DG[i] = (gamma_opt[i] - gamma_opt[mod1(i-1,Nx)])*estimate.dt/mesh.dx
end

# Diffusion
D = zero(L)
@unpack etacont_init, etacont = estimate
diffusion!(Posteriori(), gamma_opt, etacont_init, etacont, estimate.dt, mesh, D)
Dc = zero(L)
diffusion!(Posteriori(), Gc, etacont_init, etacont, estimate.dt, mesh, Dc)

pltA = []

plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topleft,
    legendfontsize=15,
    titlefontsize=21,
    guidefontsize=21,
    tickfontsize=18)

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

plot!(mesh.x, D, label="uzawa", lw=2)
plot!(mesh.x, Dexact, label="from theory", lw=2)
plot!(mesh.x, Dc, label="Gc", lw=2)
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

display(plot(pltA..., layout=(2,2), size=(1600, 1200)))

# plot(mesh.x, optsol.Gcgap)

window = 1:10
pltB=[]
plt1 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
scatter!(mesh.x[window], estimate.uinit[window], label="u0")
push!(pltB, plt1)

plt2 = plot(size=(900, 600), margin=0.5Plots.cm, legend=:topright,
legendfontsize=15,
titlefontsize=21,
guidefontsize=21,
tickfontsize=18)
plot!(mesh.x[window], optsol.gamma_opt[window], label="uzawa")
scatter!(mesh.x[window], Gexact[window], label="from theory", marker=:square)
scatter!(mesh.x[window], Gc[window], label="centred", marker=:diamond)
scatter!(mesh.x[window], Gcont[window], label="G(u)", marker=:cross)
xlabel!("x")
ylabel!("G")
title!("Entropy flux in the constant area")
push!(pltB, plt2)

display(plot(pltB..., layout=(2,1), size=(1200, 1200)))