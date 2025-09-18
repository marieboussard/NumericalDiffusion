using ForwardDiff

function f(x::Real, fvec::AbstractVector, xvec::AbstractVector)
    N = length(fvec)
    for k in 1:N
        if (xvec[k] <= x) && (x <= xvec[k+1])
            return fvec[k]
        end
    end
    fvec[N]
end


function z0(x::Real, fvec::AbstractVector, xvec::AbstractVector, alpha::AbstractVector, beta::AbstractVector)
    N = length(alpha)
    for k in 1:N
        if (xvec[k] <= x) && (x <= xvec[k+1])
            return -fvec[k]*0.5*x^2 + alpha[k]*x + beta[k]
        end
    end
    -fvec[N]*0.5*x^2 + alpha[N]*x + beta[N]
end

z0p(x::Real, fvec::AbstractVector, xvec::AbstractVector, alpha::AbstractVector, beta::AbstractVector) = ForwardDiff.derivative(t -> z0(t, fvec, xvec, alpha, beta), x)


N = 5
fvec = (2*rand(N) .- 1)*10
xvec = LinRange(0, 1, N+1)

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


xthin = LinRange(0.0, 1.0, 1000)
fthin = [f(xi, fvec, xvec) for xi in xthin]
z0vec = [z0(xi, fvec, xvec, alpha, beta) for xi in xthin]
z0pvec = [z0p(xi, fvec, xvec, alpha, beta) for xi in xthin]


using CairoMakie 
fig = Figure(size=(800, 1500))
ax = Axis(fig[1,1], xlabel="x", title="f(x)")
lines!(ax, xthin, fthin)
ax2 = Axis(fig[2,1], xlabel="x", title="z0(x)")
lines!(ax2, xthin, z0vec)
ax3 = Axis(fig[3,1], xlabel="x", title="z0'(x)")
lines!(ax3, xthin, z0pvec)
fig