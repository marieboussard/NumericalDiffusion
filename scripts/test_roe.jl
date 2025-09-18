using NumericalDiffusion
include("./optimisation/beta_slopes.jl")

Nx = 300
xmin, xmax = -2.0, 2.0
t0, tf = 0.0, 10.0
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

function u0(x::Real)
    if x < -3/2
        return one(x)
    elseif x < -1/2
        return -2*x-2.0
    elseif x < 1
        return -one(x)
    else
        return one(x)
    end
end

u0(x::AbstractVector) = -0.1*u0.(x)

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)
time_scheme = Euler()
space_scheme = Roe()

sol = solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false));
# estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()))
# @show sum(estimate.L)
# optsol = compute_entropic_G(params, equation)

# DG = zero(optsol.G_newt)
# DGc = zero(optsol.Gc)
# for j in 1:Nx
#     DG[j] = sol.dt/ mesh.dx*(optsol.G_newt[mod1(j, Nx)] - optsol.G_newt[mod1(j-1, Nx)])
#     DGc[j] = sol.dt/ mesh.dx*(optsol.Gc[mod1(j, Nx)] - optsol.Gc[mod1(j-1, Nx)])
# end


using CairoMakie
fig = Figure()
# fig = Figure(size = (800, 1500))
ax = Axis(fig[1,1], title = get_name(time_scheme)*" + "*get_name(space_scheme), xlabel="x", ylabel="u")
lines!(ax, mesh.x, sol.uinit, label="uinit", color=:navy)
lines!(ax, mesh.x, sol.u, label="t = "*string(sol.t), color=:tomato)
axislegend(ax, position = :rb)

# ax2 = Axis(fig[2,1], title="Numerical Diffusion")
# lines!(ax2, mesh.x, optsol.D_newt, label="Gopt")
# lines!(ax2, mesh.x, optsol.Dc, label="Gc")
# axislegend(ax2, position=:rb)

# ax3 = Axis(fig[3,1], title="DG")
# lines!(ax3, mesh.x, estimate.l, label="l")
# lines!(ax3, mesh.x, DG, label="DG")
# lines!(ax3, mesh.x, DGc, label="DGc")
# lines!(ax3, mesh.x, estimate.L, label="L")
# axislegend(ax3, position=:rb)

fig

# Nxlog = LinRange(1, 1.7, 4)
# #Nxlog = LinRange(1,1.5,4)
# Nxvec = Int.(floor.(10 .^ (Nxlog)))
# alpha = 1
# weights = AbsWeights(alpha)
# bound_mode = SingleBound()
# tf = 0.1

# time_scheme_vec = [Euler()]
# space_scheme_vec = [Roe()]

# fig2 = plot_beta_slopes(Nxvec, equation, time_scheme_vec, space_scheme_vec; tf=tf)