using NumericalDiffusion
include("./optimisation/beta_slopes.jl")

Nx = 50
xmin, xmax = -3.0, 3.0
t0, tf = 0.0, 0.2
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# function u0(x::Real)
#     if x < -1
#         return one(x)
#     elseif x < 1/2
#         return -one(x)
#     elseif x < 3/2
#         return 2*x - 2.0
#     else
#         return one(x)
#     end
# end

# function u0(x::Real)
#     if x < -3/2
#         return -one(x)
#     elseif x < -1/2
#         return 2*x+2.0
#         #return 3*x+7/2
#     elseif x < 1
#         return one(x)
#         #return 2*one(x)
#     else
#         return -one(x)
#     end
# end

# u0(x::AbstractVector) = -min.(u0.(x), 1)

function u0(x::Real)
    if x <= 0
        return min(0.5*x, -1)
    else
        return max(1, 0.5*x)
    end
end

u0(x::AbstractVector) = u0.(x)

# Exact solution 
function uexact(x::Real, t::Real)
    if x <= -t - 3/2
        return -1
    elseif x <= t - 1/2
        return (2*x+2)/(2*t+1)
    elseif x <= 1
        return 1
    else
        return -1
    end
end

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
time_scheme = Euler()
space_scheme = Roe()


sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false));

uevec = [uexact(xi, sol.t) for xi in sol.params.mesh.x]


# estimate = quantify_diffusion(sol, Posteriori(AsymmetricMD()))
# estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()))
# @show sum(estimate.L)
optsol = compute_entropic_G(params, equation; time_scheme=time_scheme, space_scheme=space_scheme)

DG = zero(optsol.G_newt)
DGc = zero(optsol.Gc)
for j in 1:Nx
    DG[j] = sol.dt/ mesh.dx*(optsol.G_newt[mod1(j, Nx)] - optsol.G_newt[mod1(j-1, Nx)])
    DGc[j] = sol.dt/ mesh.dx*(optsol.Gc[mod1(j, Nx)] - optsol.Gc[mod1(j-1, Nx)])
end



using CairoMakie
fig = Figure(size = (800, 1500))
ax = Axis(fig[1,1], title = get_name(time_scheme)*" + "*get_name(space_scheme), xlabel="x", ylabel="u")
lines!(ax, mesh.x, sol.uinit, label="uinit", color=:navy)
lines!(ax, mesh.x, sol.u, label="t = "*string(round(sol.t, digits=2)), color=:tomato)
#lines!(ax, mesh.x, uevec, label="exact", color=:green)
axislegend(ax, position = :rb)

# ax2 = Axis(fig[2,1], title="Numerical Diffusion")
# lines!(ax2, mesh.x, estimate.D, label="D")
# axislegend(ax2, position=:rb)

ax2 = Axis(fig[2,1], title="Numerical Diffusion")
lines!(ax2, mesh.x, optsol.D_newt, label="Dopt", color=:tomato)
lines!(ax2, mesh.x, optsol.Dc, label="Dc", color=:navy)
axislegend(ax2, position=:rb)

ax3 = Axis(fig[3,1], title="DG")
lines!(ax3, mesh.x, optsol.l, label="l", color=:black)
lines!(ax3, mesh.x, DG, label="DG", color=:tomato)
lines!(ax3, mesh.x, DGc, label="DGc", color=:navy)
lines!(ax3, mesh.x, optsol.L, label="L", color=:green)
axislegend(ax3, position=:rb)

fig

Nxlog = LinRange(1, 2.1, 10)
#Nxlog = LinRange(1,1.5,4)
Nxvec = Int.(floor.(10 .^ (Nxlog)))
alpha = 1
weights = AbsWeights(alpha)
bound_mode = SingleBound()
tf = 0.1

time_scheme_vec = [Euler()]
space_scheme_vec = [Roe()]

fig2 = plot_beta_slopes(Nxvec, equation, time_scheme_vec, space_scheme_vec; tf=tf)

fig2