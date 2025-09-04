using NumericalDiffusion
using UnPack
using LinearAlgebra

"""
    solve_poisson(; x0::Real=0, xend::Real=1; u0::Real=0, uend::Real=0)
solve Poisson equation with Dirichlet BCs u(`x0`)=`u0` and u(`xend`)=`uend`
"""
function solve_poisson!(u::AbstractVector, x::AbstractVector, dx::Real; F::AbstractVector=zero(u), x0::Real=0, xend::Real=1, u0::Real=0, uend::Real=0)

    N = length(u)

    # A = (uend - u0) / (xend - x0)
    # B = (u0 * xend - uend * x0) / (xend - x0)
    # @. u = A * x + B

    #@unpack x, dx = mesh
    A = -1 / dx^2 * zeros(eltype(u), (N, N))

    for i in 1:N
        for j in 1:N
            if i == j
                A[i, j] = -2
            elseif i == j + 1
                A[i, j] = 1
            elseif i + 1 == j
                A[i, j] = 1
            end
        end
    end

    A .*= -1 / dx^2

    #Fbc = copy(view(F, 2:N-1))
    Fbc = copy(F)
    Fbc[1] += u0 / dx^2
    Fbc[end] += uend / dx^2

    #ushort = view(u, 2:N-1)

    mul!(u, inv(A), Fbc)

    #u[1], u[N] = u0, uend

end


function create_initial_condition!(z::AbstractVector, w::AbstractVector, mesh::OneDMesh, q::AbstractVector; dt::Real=1)

    @unpack x = mesh

    N = length(q)
    qpos = q .>= 0

    w[.!qpos] .= zero(eltype(w))
    w[qpos] .= q[qpos]
    z[qpos] .= zero(eltype(z))

    # Resolving Poisson equation separately on each interval
    start_vec = Int[]
    end_vec = Int[]

    for j in 1:N
        if q[j] < 0 && q[mod1(j - 1, N)] >= 0
            push!(start_vec, j)
        end
        if q[j] < 0 && q[mod1(j + 1, N)] >= 0
            push!(end_vec, j)
        end
    end

    for k in eachindex(start_vec)
        #x0, xend = x[mod1(start_vec[k] - 1, N)], x[mod1(end_vec[k] + 1, N)]
        F = -1 / dt^2 * view(q, start_vec[k]:end_vec[k])
        solve_poisson!(view(z, start_vec[k]:end_vec[k]), view(x, start_vec[k]:end_vec[k]), mesh.dx; F=F)#; x0=x0, xend=xend)
    end

end


Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

time_scheme = Euler()
space_scheme = Rusanov()

bound_mode = SingleBound()

weights = AbsWeights(0)


# Finite volumes resolution
sol = NumericalDiffusion.solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false))

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()))
@unpack uinit, u, l, L = estimate

Gc, A, b, W = init_optim_components(bound_mode, estimate, weights)

q = b - A * Gc
w0 = zero(b)
z0 = zero(b)


@time create_initial_condition!(z0, w0, mesh, q; dt=estimate.dt)

# Run newton with it 
M = A * inv(W) * A'
@time pend, wend, niter = newton_lcp(M, -q, z0, w0; maxiter=1000)
# We get the flux back 
Gend = Gc - inv(W) * A' * pend
newton_residual = norm(max.(0.0, A * Gend .- b))
@show newton_residual
@unpack etacont_init, etacont = estimate
Dend = zero(L)
diffusion!(Posteriori(), Gend, etacont_init, etacont, estimate.dt, params.mesh, Dend)

# Compare to the final result
@time refsol = compute_entropic_G(params, equation; bound_mode=bound_mode, weights=weights)
zref = refsol.p_newt
wref = refsol.b - refsol.A * refsol.G_newt

# Plot

using CairoMakie

fig = Figure(resolution=(1000, 1000))
ax = Axis(fig[1, 1], title="z", xlabel="x")# =get_name(time_scheme) * " + " * get_name(space_scheme), xlabel="x")

lines!(ax, mesh.x, z0, color=:tomato, label="z0")
scatter!(ax, mesh.x, z0, color=:tomato)
lines!(ax, mesh.x, zref, color=:navy, label="zref")
scatter!(ax, mesh.x, zref, color=:navy)
#lines!(ax, mesh.x, pend, color=:green, label="zend")
scatter!(ax, mesh.x, pend, color=:green, label="zend")
axislegend(position=:ct)

ax2 = Axis(fig[2, 1], title="w", xlabel="x")
lines!(ax2, mesh.x, w0, color=:tomato, label="w0")
scatter!(ax2, mesh.x, w0, color=:tomato)
lines!(ax2, mesh.x, wref, color=:navy, label="wref")
scatter!(ax2, mesh.x, wref, color=:navy)
#lines!(ax2, mesh.x, wend, color=:green, label="wend")
scatter!(ax2, mesh.x, wend, color=:green, label="wend")
axislegend(position=:cb)

fig

fig2 = Figure()
ax = Axis(fig2[1, 1], title="Numerical Diffusion", xlabel="x")
lines!(ax, mesh.x, Dend, color=:green, label="Dend")
scatter!(ax, mesh.x, Dend, color=:green)
lines!(ax, mesh.x, refsol.D_newt, color=:navy, label="Dref")
scatter!(ax, mesh.x, refsol.D_newt, color=:navy)
axislegend(position=:lb)

fig2

fig3 = Figure()
ax = Axis(fig3[1, 1], title="Numerical entropy flux", xlabel="x")
lines!(ax, mesh.x, refsol.G_newt, color=:navy, label="Gref")
scatter!(ax, mesh.x, Gend, color=:green, label="Gend")
#scatter!(ax, mesh.x, Gend, color=:green)
#lines!(ax, mesh.x, refsol.G_newt, color=:navy, label="Gref")

axislegend(position=:lb)

fig3