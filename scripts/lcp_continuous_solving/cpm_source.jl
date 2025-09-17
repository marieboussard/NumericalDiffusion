using ForwardDiff
using Integrals
using CairoMakie

# beta = zeros(N)
# beta[1] = 10
# for k in 2:N
#     beta[k] = beta[k-1] +0.5*x[k]^2*(f[k-1]-f[k])
# end

# alpha = zeros(N)
# alpha[N] = f[N]/2 + beta[1] - beta[N]
# for k in N-1:-1:1
#     alpha[k] = alpha[k+1] - (f[k+1]-f[k])*x[k+1]
# end

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

z(x::Real, fvec::AbstractVector, xvec::AbstractVector, alpha::AbstractVector, beta::AbstractVector) = max(zero(x), z0(x, fvec, xvec, alpha, beta))

zp(x::Real, fvec::AbstractVector, xvec::AbstractVector, alpha::AbstractVector, beta::AbstractVector) = ForwardDiff.derivative(t -> z(t, fvec, xvec, alpha, beta), x)

g(x::Real, fvec::AbstractVector, xvec::AbstractVector, alpha::AbstractVector, beta::AbstractVector) = 0.5 * (zp(x, fvec, xvec, alpha, beta))^2 - f(x, fvec, xvec)*z(x, fvec, xvec, alpha, beta)


function J(y::Real, fvec::AbstractVector, xvec::AbstractVector; plot=false)
    N = length(fvec)

    # Compute coefficients of z0 on each interval
    beta = zeros(N)
    beta[1] = y
    for k in 2:N
        beta[k] = beta[k-1] +0.5*xvec[k]^2*(fvec[k-1]-fvec[k])
    end
    alpha = zeros(N)
    alpha[N] = fvec[N]/2 + beta[1] - beta[N]
    for k in N-1:-1:1
        alpha[k] = alpha[k+1] - (fvec[k+1]-fvec[k])*xvec[k+1]
    end

    domain = (0.0, 1.0)
    prob = IntegralProblem((x, p)->g(x, fvec, xvec, alpha, beta), domain)
    sol = solve(prob, HCubatureJL(); reltol = 1e-3, abstol = 1e-3)

    if plot

        xvec_thin = LinRange(0, 1, 100)
        zvec = [z(x, fvec, xvec, alpha, beta) for x in xvec_thin]

        temp_fig = Figure()
        ax = Axis(temp_fig[1,1], xlabel="x", ylabel="z(x)", title="beta = "*string(y))
        lines!(ax, xvec_thin, zvec)

        return temp_fig, sol.u, alpha, beta

    else

        return sol.u, alpha, beta

    end

end


N = 5
fvec = (2*rand(N) .- 1)*10
@show sum(fvec)
#fvec = [-1.0]
xvec = LinRange(0, 1, N+1)

yvec = LinRange(-1.0, 1, 20)
Jvec = zero(yvec)
figvec = Figure[]
for k in eachindex(yvec)
    t_fig, Jvec[k], alpha, beta = J(yvec[k], fvec, xvec; plot=true)
    push!(figvec, t_fig)
end

#Jvec = [J(y, fvec, xvec)[1] for y in yvec]

# Treshold 
x_thin = LinRange(0, 1, 1000)
Jbar, alphabar, betabar = J(0.0, fvec, x_thin)
beta_t = -Inf
for j in eachindex(x_thin)
    if f(x_thin[j], fvec, xvec) >= 0
        test = -z0(x_thin[j], fvec, xvec, alphabar, betabar)
        if test > beta_t
            global beta_t = test
        end
    end
end

@show beta_t

z0bar_vec = [z0(xi, fvec, xvec, alphabar, betabar) for xi in x_thin]
fvec_thin = [f(xi, fvec, xvec) for xi in x_thin]
fpos = fvec_thin .>= 0

fig = Figure(size=(800,1000))
ax = Axis(fig[1,1], xlabel="x", title="f(x)")
lines!(ax, x_thin, fvec_thin)
scatter!(ax, x_thin[fpos], fvec_thin[fpos], color=:tomato)

axb = Axis(fig[2,1], xlabel="x", title="z0(x,beta=0)")
lines!(axb, x_thin, z0bar_vec)
scatter!(axb, x_thin[fpos], z0bar_vec[fpos], color=:tomato)
lines!(axb, [0.0, 1.0], [-beta_t, -beta_t], color=:green)




fig2 = Figure()

ax2 = Axis(fig2[1,1], xlabel="beta", ylabel="J(beta)")
lines!(ax2, yvec, Jvec)
scatter!(ax2, yvec, Jvec)
lines!(ax2, [beta_t, beta_t], [minimum(Jvec), maximum(Jvec)])
fig2



# z0p(t) = ForwardDiff.derivative(z0, t)
# z0pp(t) = ForwardDiff.derivative(z0p, t)

# xvec = LinRange(0,1, 100)
# z0vec = z0.(xvec)#, alpha, beta)
# z0pvec = z0p.(xvec)
# z0ppvec = z0pp.(xvec)


# using CairoMakie 
# fig = Figure(size=(800, 1500))

# ax1 = Axis(fig[1,1], xlabel="x", ylabel="z0")
# scatter!(ax1, xvec, z0vec)

# ax2 = Axis(fig[2,1], xlabel="x", ylabel="z0'")
# scatter!(ax2, xvec, z0pvec)

# ax2 = Axis(fig[3,1], xlabel="x", ylabel="z0''")
# scatter!(ax2, xvec, z0ppvec)

# fig