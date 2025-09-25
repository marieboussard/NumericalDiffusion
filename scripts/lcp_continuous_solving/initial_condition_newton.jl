using NumericalDiffusion
using UnPack
include("exact_poisson.jl")

function solve_poisson_dirichlet_exact!(u::AbstractVector, fvec::AbstractVector, xvec::AbstractVector)
    N = length(fvec)

    # Compute coefficients of the piecewise quadratic solution
    alpha = zeros(N)
    beta = zeros(N)
    alpha[1] = 0.5*(fvec[1]*xvec[1]^2-fvec[N]*xvec[N+1]^2)
    for j in 1:N-1
        alpha[1] += (fvec[j+1]-fvec[j])*xvec[j+1]*(xvec[N+1]-0.5*xvec[j+1])
    end
    alpha[1] *= 1.0/(xvec[1]-xvec[N+1])
    for j in 1:N-1
        alpha[j+1] = alpha[j] + (fvec[j+1]-fvec[j])*xvec[j+1]
    end
    beta[1] = 0.5*fvec[1]*xvec[1]^2 - alpha[1]*xvec[1]
    for j in 1:N-1
        beta[j+1] = beta[j] - 0.5*xvec[j+1]^2*(fvec[j+1]-fvec[j])
    end

    # Fill uj ≈ u(xj)
    for j in 1:N
        xj = 0.5*(xvec[j]+xvec[j+1])
        u[j] = -fvec[j]*0.5*xj^2 + alpha[j]*xj + beta[j]
    end

end

function create_initial_condition!(zvec::AbstractVector, wvec::AbstractVector, mesh::OneDMesh, q::AbstractVector; dt::Real=1, CFL_factor::Real=0.5)

    @unpack interfaces = mesh

    N = length(q)
    qpos = q .>= 0

    wvec[.!qpos] .= zero(eltype(wvec))
    wvec[qpos] .= q[qpos]
    zvec[qpos] .= zero(eltype(zvec))

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

        F = -view(q, start_vec[k]:end_vec[k])/dt^2

        xshort = view(interfaces, start_vec[k]:end_vec[k]+1)
        zshort = view(zvec, start_vec[k]:end_vec[k])

        solve_poisson_dirichlet_exact!(zshort, F, xshort)


    end

end

function create_initial_condition_max!(zvec::AbstractVector, wvec::AbstractVector, mesh::OneDMesh, q::AbstractVector; dt::Real=1, CFL_factor::Real=0.5, betaone=0.0)

    fvec = -q/(dt^2)
    xvec = mesh.interfaces
    source = CPMSource(fvec, xvec)
    problem = PoissonProblem(source, Periodic())
    solve_poisson_cpm!(zvec, wvec, mesh.x, problem; betaone=betaone)
    wvec .*= dt^2
    @show problem.source.alpha
    @show problem.source.beta

end

Nx = 200
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)
equation = BurgersArticle

time_scheme = Euler()
space_scheme = Rusanov()
bound_mode = SingleBound()
weights = AbsWeights(1)

# Finite volumes resolution
sol = NumericalDiffusion.solve(equation, params, time_scheme, space_scheme; log_config=LogConfig(true, false, true, false, false))

# Multidimensional bounds for ΔG
estimate = quantify_diffusion(sol, PrioriMultidim(AsymmetricMD()))
@unpack uinit, u, l, L = estimate

Gc, A, b, W = init_optim_components(bound_mode, estimate, weights)

q = b - A * Gc

# FIRST METHOD WITH HOLES
println("COMPUTING INITIAL CONDITION WITH HOLES")
winit = zero(b)
zinit = zero(b)
create_initial_condition!(zinit, winit, mesh, q; dt=estimate.dt, CFL_factor=params.CFL_factor)

# SECOND METHOD WITH MAX
println("COMPUTING INITIAL CONDITION WITH MAX")
wmax = zero(b)
zmax = zero(b)
betaone = 5.0
create_initial_condition_max!(zmax, wmax, mesh, q; dt=estimate.dt, CFL_factor=params.CFL_factor, betaone=betaone)

# REFERENCE METHOD
println()
println("REFERENCE UZAWA + NEWTON")
refsol = compute_entropic_G(params, equation; bound_mode=bound_mode, weights=weights)
zref = refsol.p_newt
wref = refsol.b - refsol.A * refsol.G_newt

using CairoMakie
fig = Figure(size=(1000, 1000))

ax = Axis(fig[1, 1], title="z(x)", xlabel="x")
lines!(ax, mesh.x, zinit, color=:tomato, label="holes")
scatter!(ax, mesh.x, zinit, color=:tomato)
lines!(ax, mesh.x, zmax, color=:navy, label="max")
scatter!(ax, mesh.x, zmax, color=:navy)
lines!(ax, mesh.x, zref, color=:green, label="zref")
scatter!(ax, mesh.x, zref, color=:green)
axislegend(ax, position=:lt)

ax2 = Axis(fig[2, 1], title="w", xlabel="x")
lines!(ax2, mesh.x, winit, color=:tomato, label="holes")
scatter!(ax2, mesh.x, winit, color=:tomato)
lines!(ax2, mesh.x, wmax, color=:navy, label="max")
scatter!(ax2, mesh.x, wmax, color=:navy)
lines!(ax2, mesh.x, wref, color=:green, label="wref")
scatter!(ax2, mesh.x, wref, color=:green)
axislegend(ax2, position=:lt)

fig