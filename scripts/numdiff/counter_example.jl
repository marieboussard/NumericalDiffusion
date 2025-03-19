using BenchmarkTools
include("../../src/numdiff/include_file.jl")

# PARAMETRIZED INITIAL DATA

const a1 = 1.0
const b1 = 2.0
const a2 = 2.0
const b2 = 4.5

xmin, xmax = -b1/a1, b2/a2
M = 100
Nx = Int((1 + floor(M/(a1*b2+a2*b1)))*(a2*b1+a1*b2))
t0, tf = 0.0, 0.2
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# INITIAL CONDITION
u0(x::Real) = x<=0 ? -a1*x-b1 : -a2*x+b2
u0(x::AbstractArray) = u0.(x)

# EXACT SOLUTION
function uexact_fun(x::Real, t::Real)
    if t >= min(1/a1, 1/a2)
        @warn "Warning: This solution is not valid for t ≥ min(1/a1, 1/a2)"
    end
    if x <= -b1*t
        return (a1*x+b1) / (a1*t-1)
    elseif x >= b2*t
        return (a2*x-b2) / (a2*t-1)
    else
        return x/t
    end
end
uexact = [uexact_fun(xi, tf) for xi in mesh.x]

# FINITE VOLUMES RESOLUTION
equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
sol = FiniteVolumes.solve(equation, params, Euler(), Roe(); maxiter=1, log_config=LogConfig(true, false, true, false));

using Plots
plot(mesh.x, sol.uinit, label="uinit")
plot!(mesh.x, uexact, label="exact")
display(plot!(mesh.x, sol.u, label="t = "*string(sol.t)))

# QUANTIFICATION OF NUMERICAL DIFFUSION
estimate_priori = quantify_diffusion(sol, PrioriMultidim(MaxMD()); name="priori multidim");
estimate_posteriori = quantify_diffusion(sol, Posteriori(AsymmetricMD()); name="asymmetric posteriori");
estimate_prioristd = quantify_diffusion(sol, Priori(AsymmetricMD()); name="priori standard");

# CHECKING COUNTER EXAMPLE CONDITIONS

println("min(M-m):")
println(minimum(estimate_posteriori.M .- estimate_posteriori.m))
println("min(L-l)")
println(minimum(estimate_priori.L .- estimate_priori.l))
println("sum(L)")
println(sum(estimate_priori.L))
println("sum(l)")
println(sum(estimate_priori.l))

plot(mesh.x, estimate_priori.l, label="l")
display(plot!(mesh.x, estimate_priori.L, label="L"))

plot(mesh.x, estimate_posteriori.m, label ="m")
display(plot!(mesh.x, estimate_posteriori.M, label ="M"))

plot(mesh.x, estimate_posteriori.D, label="posteriori")
plot!(mesh.x, estimate_priori.D, label="priori multidim")
plot!(mesh.x, estimate_prioristd.D, label = "priori standard")
plot!(mesh.x, estimate_priori.Dlow, label="Dlow")
display(title!("Numerical Diffusion"))