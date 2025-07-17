include("../../src/optimisation/uzawa_newton.jl")
using Plots

Nx = 50
xmin = -2.0
xmax = 2.0
t0 = 0.0
tf = 0.5
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

function u0(x::Real)
    if x<-1
        return 1.0
    elseif x <1
        return -1.0
    else
        return 1.0
    end
end
u0(x::AbstractVector) = u0.(x)

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)

# Finite volumes resolution
sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, false, true, false, false));

plot(params.mesh.x, sol.uinit)
display(plot!(params.mesh.x, sol.u))

estimate = quantify_diffusion(sol, Posteriori(AsymmetricMD()));

@show sum(estimate.D)