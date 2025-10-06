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


function solve_with_matrix(fvec::AbstractVector, xvec::AbstractVector)

    M = length(fvec)
    A = zeros(2*M,2*M)
    F = zeros(2*M)

    # Dirichlet condition on left boundary 
    A[1,1] = xvec[1]
    A[1,2] = 1.0
    F[1] = 0.5*fvec[1]*xvec[1]^2

    for i in 1:M-1

        # Continuity of z
        A[2*i, 2*i-1] = xvec[i+1]
        A[2*i, 2*i] = 1.0
        A[2*i, 2*i+1] = -xvec[i+1]
        A[2*i, 2*i+2] = -1.0
        F[2*i] = 0.5*(fvec[i]+fvec[i+1])*xvec[i+1]^2

        # Continuity of z'
        A[2*i+1, 2*i-1] = 1.0
        A[2*i+1, 2*i+1] = -1.0
        F[2*i+1] = (fvec[i]-fvec[i+1])*xvec[i+1]

    end

    # Dirichlet condition on right boundary
    A[2*M, 2*M-1] = xvec[M+1]
    A[2*M, 2*M] = 1.0
    F[2*M] = 0.5*fvec[M]*xvec[M+1]^2

    inv(A)*F

end

function solve_block(fvec::AbstractVector, xvec::AbstractVector)

    N = length(fvec)
    fpos = fvec .>=0

    # Resolving Poisson equation separately on each interval
    start_vec = Int[]
    end_vec = Int[]

    for j in 1:N
        if fvec[j] >= 0 && (j == 1 || fvec[j-1] < 0)
            push!(start_vec, j)
        end
        if fvec[j] >= 0 && (j == N || fvec[j+1] < 0)
            push!(end_vec, j)
        end
    end

    @show fvec
    @show start_vec, end_vec

    coeffs = []


    for k in eachindex(start_vec)

        fshort = view(fvec, start_vec[k]:end_vec[k])

        xshort = view(xvec, start_vec[k]:end_vec[k]+1)
        #zshort = view(zvec, start_vec[k]:end_vec[k])

        push!(coeffs, solve_with_matrix(fshort, xshort))

    end

    return start_vec, end_vec, coeffs

end

function cont_sol(x::Real, fvec::AbstractVector, xvec::AbstractVector, start_vec::AbstractVector, end_vec::AbstractVector, coeffs::AbstractVector)
    K = length(start_vec)

    for k in 1:K
        if xvec[start_vec[k]]<= x && xvec[end_vec[k]+1] >= x
            for i in start_vec[k]:end_vec[k]
                if xvec[i] <= x && xvec[i+1] >= x
                    return -0.5*fvec[i]*x^2 + coeffs[k][2*i-1]*x + coeffs[k][2*i], 0.0
                end
            end
        end

    end

    N = length(fvec)
    for j in 1:N
        if xvec[j] <= x && xvec[j+1] >= x 
            return 0.0, -fvec[j]*Delta^2#*(0.5*(xvec[end]-xvec[begin])/N)^2
        end 
    end

    return 0.0, 0.0

end


Delta = 1.0



Nxvec = [20, 50, 100, 200]
xmin, xmax = -1, 1
t0, tf = 0.0, 0.4
#dtN = 0.5*(xmax-xmin)/N

# # Random source term
# function gen_random_q!(q::Vector{Float64}, N::Int)
#     q .= (2 * rand(N) .- 1.0) .* 1
# end
# N = 5
# qN = zeros(N)
# gen_random_q!(qN, N)
# while sum(qN) < 0
#     gen_random_q!(qN, N)
# end

N = 2
qN = [-1.0, 2.0]


# 2 # Continuous solution
Nx = 100

mesh = OneDMesh(Nx, xmin, xmax)
dt = 0.5 * mesh.dx

#q = [qN[Int(floor(i / Nx * N) + 1)] for i in 0:Nx-1]
fvec = -qN/Delta^2#/(0.5*(xmax-xmin)/N)^2
xvec = LinRange(xmin, xmax, N+1)
#solve_continuous_cp_dirichlet!(zcont, wcont, mesh, q; dt=dt)
start_vec, end_vec, coeffs = solve_block(fvec,xvec)


wcont = zeros(Nx)
zcont = zeros(Nx)
for i in 1:Nx
    zcont[i], wcont[i] = cont_sol(mesh.x[i], fvec, xvec, start_vec, end_vec, coeffs)
end




# Visualization
using CairoMakie
fig = Figure(size=(1000, 1000))

ax = Axis(fig[1, 1], xlabel="x", title="z")

#lines!(ax, mesh.x, zcont, color=:navy, label="continuous")
#scatter!(ax, mesh.x, zcont, color=:navy, label="continuous")

ax2 = Axis(fig[2, 1], xlabel="x", title="w")
#lines!(ax2, mesh.x, wcont, color=:navy, label="continuous")



function iterate_discret(ax, ax2, Nxvec::AbstractVector, qN, xmin=-1, xmax=1, t0=0.0, tf=0.4)

    Nvec = length(Nxvec)


    for k in 1:Nvec

        Nx = Nxvec[k]

        mesh = OneDMesh(Nx, xmin, xmax)
        dt = 0.5 * mesh.dx


        wcont = zeros(Nx)
        zcont = zeros(Nx)
        for i in 1:Nx
            zcont[i], wcont[i] = cont_sol(mesh.x[i], fvec, xvec, start_vec, end_vec, coeffs)
        end

        

        lines!(ax, mesh.x, zcont, label="Nx="*string(Nx))
        #scatter!(ax, mesh.x, zcont, color=:navy, label="continuous")

        lines!(ax2, mesh.x, wcont, label="continuous")



        q = -1*[qN[Int(floor(i / Nx * N) + 1)] for i in 0:Nx-1]

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
        #M .*= (dt / mesh.dx)^2
        M .*= (Delta / mesh.dx)^2

        @show maximum(abs.(M*zcont .- wcont .+ q))

        # 1 # Solution of related LCP problem
        # z0vec = zeros(Nx)
        # w0vec = zeros(Nx)
        z0vec = zcont 
        w0vec = wcont
        zvec, wvec, niter = newton_lcp(M, q, z0vec, w0vec; printing=true, maxiter=1000)
        # result = solve!(LCP(M, q); max_iter=100000)
        # @show result
        # zvec = result.sol
        # wvec = q .+ M * zvec


        scatter!(ax, mesh.x, zvec, label="Nx="*string(Nx))
        scatter!(ax2, mesh.x, wvec, label="Nx="*string(Nx))
    

    end
end

iterate_discret(ax, ax2, Nxvec, qN, xmin, xmax, t0, tf)

axislegend(ax, position=:ct)
axislegend(ax2, position=:lt)
fig