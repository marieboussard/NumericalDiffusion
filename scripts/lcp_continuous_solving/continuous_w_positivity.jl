using ForwardDiff

# # f = constant
# z0(x, beta, f) = -f*x^2 / 2 + f / 2 * x + beta
# z(x, beta, f) = max(z0(x, beta, f), 0.0)
# zp(x, beta, f) = ForwardDiff.derivative(t -> z(t, beta, f), x)
# zpp(x, beta, f) = ForwardDiff.derivative(t -> zp(t, beta, f), x)

# N = 20
# f = 1
# beta = 1
# xvec = LinRange(0, 1, N)
# z0vec = z0.(xvec, beta, f)
# zvec = z.(xvec, beta, f)
# zpvec = zp.(xvec, beta, f)
# zppvec = zpp.(xvec, beta, f)
# wvec = -zppvec .- f

# using CairoMakie
# fig = Figure()
# ax = Axis(fig[1, 1], title="f="*string(f)*", beta="*string(beta), xlabel="x")

# lines!(ax, xvec, z0vec, color=:green, label="z0")
# scatter!(ax, xvec, z0vec, color=:green)

# lines!(ax, xvec, wvec, color=:tomato, label="w")
# scatter!(ax, xvec, wvec, color=:tomato)

# lines!(ax, xvec, zvec, label="z", color=:navy)
# scatter!(ax, xvec, zvec, color=:navy)

# # lines!(ax, xvec, zpvec, label="z'")
# # scatter!(ax, xvec, zpvec)

# # lines!(ax, xvec, zppvec, label="z''")
# # scatter!(ax, xvec, zppvec)

# axislegend(ax, position=:lc)
# fig


# f = sinus(2pix)

f(x, k) = sin(2*pi*k*x)

z0(x, beta, k) = sin(2*pi*k*x)/(2*pi*k)^2 + beta
z(x, beta, k) = max(z0(x, beta, k), 0.0)
zp(x, beta, k) = ForwardDiff.derivative(t -> z(t, beta, k), x)
zpp(x, beta, k) = ForwardDiff.derivative(t -> zp(t, beta, k), x)

N = 20
beta = 0
k = 1
xvec = LinRange(0, 1, N)
zvec = z.(xvec, beta, k)
zpvec = zp.(xvec, beta, k)
zppvec = zpp.(xvec, beta, k)
wvec = -zppvec .- f.(xvec, k)

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], title="f=1", xlabel="x")
lines!(ax, xvec, wvec, color=:tomato, label="w")
scatter!(ax, xvec, wvec, color=:tomato)

lines!(ax, xvec, zvec, label="z", color=:navy)
scatter!(ax, xvec, zvec, color=:navy)

axislegend(ax, position=:lc)
fig