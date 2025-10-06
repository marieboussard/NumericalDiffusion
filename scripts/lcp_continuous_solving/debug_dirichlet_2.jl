using LinearAlgebra
using LCPsolve
using CairoMakie


function solve_with_matrix(fvec::AbstractVector, xvec::AbstractVector)

    M = length(fvec) # Number of interval of the CPM function
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
        F[2*i] = 0.5*(fvec[i]-fvec[i+1])*xvec[i+1]^2

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

function cont_sol(x::Real, fvec::AbstractVector, xvec::AbstractVector, coeffs::AbstractVector)

    N = length(fvec)
    z = zero(x)

    for j in 1:N 
        if xvec[j] <= x && x < xvec[j+1]
            z = -0.5*fvec[j]*x^2 + coeffs[2*j-1]*x + coeffs[2*j]
            #return z, zero(z)
            if z >= 0.0
                return z, zero(z)
            else
                return zero(z), -fvec[j]
            end
        end
    end
    z = -0.5*fvec[N]*x^2 + coeffs[2*N-1]*x + coeffs[2*N]
    if z >= 0.0
        return z, zero(z)
    else
        return zero(z), -fvec[j]
    end
    #return z, zero(z)
end

function iterate_discret!(ax::Axis, ax2::Axis, Nxvec::AbstractVector, fvec::AbstractVector, xvec::AbstractVector, coeffs::AbstractVector, xmin=-1, xmax=1)

	N = length(fvec)
    Nvec = length(Nxvec)

    for k in 1:Nvec
        
        Nx = Nxvec[k]
        finterp = [fvec[Int(floor(i / Nx * N) + 1)] for i in 0:Nx-1]

		# Mz - w = f

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
		dx = (xmax-xmin)/Nx
		x = LinRange(xmin, xmax, Nx)
        M .*= (1.0 / dx)^2

        #@show maximum(abs.(M*zcont .- wcont .- finterp))

        # 1 # Solution of related LCP problem
        z0vec = zeros(Nx)
        w0vec = zeros(Nx)
        for i in 1:Nx
            z0vec[i], w0vec[i] = cont_sol(x[i], fvec, xvec, coeffs)
        end

        @show maximum(abs.(M*z0vec .- w0vec .- finterp))

        # z0vec = zcont 
        # w0vec = wcont
        # zvec, wvec, niter = newton_lcp(M, finterp, z0vec, w0vec; printing=true, maxiter=1000)
        result = solve!(LCP(M, -finterp), z0vec; max_iter=100000)
        @show result
        zvec = result.sol
        wvec = M * zvec .- finterp

        scatter!(ax, x, zvec, label="Nx="*string(Nx))
        scatter!(ax2, x, wvec, label="Nx="*string(Nx))
    

    end

end


# 1 # Create a random piecewise continuous data
N = 5
fvec = rand(N)*2 .- 1.0
#fvec = [-1,-2,1,2]
xmin, xmax = -1, 1
xvec = LinRange(xmin, xmax, N+1)

# 2 # CONTINUOUS OBSTACLE PROBLEM

# Compute the polynomial coefficients of associated obstacle solution
coeffs = solve_with_matrix(fvec, xvec)

# Sample this solution on a grid
Nx = 1000
x = LinRange(xmin, xmax, Nx)
wcont = zeros(Nx)
zcont = zeros(Nx)
for i in 1:Nx
    zcont[i], wcont[i] = cont_sol(x[i], fvec, xvec, coeffs)
end

# 3 # DISCRETE RESOLUTION FOR DIFFERENT RESOLUTIONS

Nxvec = [20, 50, 100]
#Nxvec = [200]

# 4 # VISUALISATION

fig = Figure(size=(1000,1000))
ax=Axis(fig[1,1], xlabel="x", title="z")
ax2=Axis(fig[2,1], xlabel="x", title="w")
lines!(ax, x, zcont, label="continuous")
lines!(ax2,x, wcont, label="continuous")

iterate_discret!(ax, ax2, Nxvec, fvec, xvec, coeffs, xmin, xmax)

axislegend(ax, position=:lt)
axislegend(ax2, position=:lt)
fig