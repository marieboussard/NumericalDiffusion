using NumericalDiffusion
using LinearAlgebra
using LCPsolve
include("initial_condition_newton.jl")


function solve_continuous_cp_dirichlet!(zvec::AbstractVector, wvec::AbstractVector, mesh::OneDMesh, q::AbstractVector; dt::Real=1)

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
        if q[j] < 0 && (j == 1 || q[j-1] >= 0)
            push!(start_vec, j)
        end
        if q[j] < 0 && (j == N || q[j+1] >= 0)
            push!(end_vec, j)
        end
    end

    for k in eachindex(start_vec)

        F = -view(q, start_vec[k]:end_vec[k]) / dt^2

        xshort = view(interfaces, start_vec[k]:end_vec[k]+1)
        zshort = view(zvec, start_vec[k]:end_vec[k])

        solve_poisson_dirichlet_exact!(zshort, F, xshort)

    end
end


Nx = 300
xmin, xmax = -1, 1
t0, tf = 0.0, 0.4
mesh = OneDMesh(Nx, xmin, xmax)
dt = 0.5 * mesh.dx

# Random source term
function gen_random_q!(q::Vector{Float64}, N::Int)
    q .= (2 * rand(N) .- 1.0) .* 1
end
N = 5
qN = zeros(N)
gen_random_q!(qN, N)
while sum(qN) < 0
    gen_random_q!(qN, N)
end


q = [qN[Int(floor(i / Nx * N) + 1)] for i in 0:Nx-1]

# Dirichlet Laplacian matrix
M = zeros(Nx, Nx)
for j in 1:Nx
    for i in 1:Nx
        if i == j
            M[i, j] = 2
        elseif (i == j + 1) || (j == i + 1)
            M[i, j] = -1
        end
    end
end
M .*= (dt / mesh.dx)^2

# 1 # Solution of related LCP problem
z0vec = zeros(Nx)
w0vec = zeros(Nx)
#zvec, wvec, niter = newton_lcp(M, q, z0vec, w0vec; printing=true, maxiter=10000)
result = solve!(LCP(M, q); max_iter=10000)
@show result
zvec = result.sol
wvec = q .+ M * zvec

# 2 # Continuous solution
wcont = zeros(Nx)
zcont = zeros(Nx)
solve_continuous_cp_dirichlet!(zcont, wcont, mesh, q; dt=dt)

# Visualization
using CairoMakie
fig = Figure(size=(1000, 1000))

ax = Axis(fig[1, 1], xlabel="x", title="z")
lines!(ax, mesh.x, zvec, color=:tomato, label="discret")
scatter!(ax, mesh.x, zvec, color=:tomato, label="discret")
lines!(ax, mesh.x, zcont, color=:navy, label="continuous")
scatter!(ax, mesh.x, zcont, color=:navy, label="continuous")
axislegend(ax, position=:ct)

ax2 = Axis(fig[2, 1], xlabel="x", title="w")
lines!(ax2, mesh.x, wvec, color=:tomato, label="discret")
scatter!(ax2, mesh.x, wvec, color=:tomato, label="discret")
lines!(ax2, mesh.x, wcont, color=:navy, label="continuous")
scatter!(ax2, mesh.x, wcont, color=:navy, label="continuous")
axislegend(ax2, position=:lt)

fig