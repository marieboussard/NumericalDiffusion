using ForwardDiff
using Integrals
using CairoMakie 


f(x::Real, kvec::AbstractVector, fvec::AbstractVector) = sum(fvec .* sin.(2*pi*kvec*x))

function z0k(x::Real, k::Real, beta::Real)
    if k==0
        return beta
    else
        return sin(2*pi*k*x)/(2*pi*k)^2 + beta
    end
end
#z0(x::Real, kvec::AbstractVector, fvec::AbstractVector, beta::Real) = sum(fvec .* z0k.(x, kvec, beta))

function z0(x::Real, kvec::AbstractVector, fvec::AbstractVector, beta::Real)
    res = zero(x)
    for j in eachindex(kvec)
        res += fvec[j]*z0k(x, kvec[j], beta)
    end
    res
end

z(x::Real, kvec::AbstractVector, fvec::AbstractVector, beta::Real) = max(0.0, z0(x, kvec, fvec, beta))
zp(x::Real, kvec::AbstractVector, fvec::AbstractVector, beta::Real) = ForwardDiff.derivative(t -> z(t, kvec, fvec, beta), x)

g(x::Real, kvec::AbstractVector, fvec::AbstractVector, beta::Real) = 0.5 * (zp(x, kvec, fvec, beta))^2 - f(x, fvec, kvec)*z(x, kvec, fvec, beta)

function J(beta::Real, kvec::AbstractVector, fvec::AbstractVector)
    domain = (0.0, 1.0)
    prob = IntegralProblem((x, p)->g(x, kvec, fvec, beta), domain)
    sol = solve(prob, HCubatureJL(); reltol = 1e-3, abstol = 1e-3)

    sol.u
end


N = 6
fvec = (2*rand(N) .- 1)*10
kvec = [i for i in 0:N-1]
betavec = LinRange(-1, 1, 30)
xvec = LinRange(0, 1, 1000)


figvec = Figure[]
for j in eachindex(betavec)
    beta = betavec[j]
    zvec = [z0(xi, kvec, fvec, beta) for xi in xvec]
    temp_fig = Figure()
    temp_ax = Axis(temp_fig[1,1], xlabel="x", ylabel="z(x)", title="beta="*string(beta))
    lines!(temp_ax, xvec, zvec, color=:tomato)
    push!(figvec, temp_fig)
end

Jvec = [J(betavec[j], kvec, fvec) for j in eachindex(betavec)]


# Treshold 
x_thin = LinRange(0, 1, 1000)
#Jbar = J(0.0, kvec, fvec)
maxi = 0.0
beta_t = -Inf
for j in eachindex(x_thin)
    if f(x_thin[j], kvec, fvec) >= 0
        test = -z0(x_thin[j], kvec, fvec, 0.0)
        if test > beta_t
            #global maxi = test
            global beta_t = test
        end
    end
end

@show beta_t

z0bar_vec = [z0(xi, kvec, fvec, 0.0) for xi in xvec]


fvec = [f(xi, kvec, fvec) for xi in xvec]
fpos = fvec .>= 0
fig = Figure(size=(800,1000))
ax = Axis(fig[1,1], xlabel="x", title="f(x)")
lines!(ax, xvec, fvec)
scatter!(ax, xvec[fpos], fvec[fpos], color=:tomato)

axb = Axis(fig[2,1], xlabel="x", title="z0(x,beta=0)")
lines!(axb, xvec, z0bar_vec)
scatter!(axb, xvec[fpos], z0bar_vec[fpos], color=:tomato)
lines!(axb, [0.0, 1.0], [-beta_t, -beta_t], color=:green)

fig

fig2 = Figure()
ax2 = Axis(fig2[1,1], xlabel="beta", ylabel="J")
lines!(ax2, betavec, Jvec, color=:tomato)
scatter!(ax2, betavec, Jvec, color=:tomato)
lines!(ax2, [beta_t, beta_t], [minimum(Jvec), maximum(Jvec)])
fig2

# fig3 = Figure()
# ax3 = Axis(fig3[1,1], xlabel="x", ylabel="z")
# lines!(ax3, xvec, zvec)
# fig3