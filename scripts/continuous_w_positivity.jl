using ForwardDiff

z0(x, beta) = x^2 / 2 + (beta - 1 / 2) * x + beta
z(x, beta) = max(z0(x, beta), 0.0)
zp(x, beta) = ForwardDiff.derivative(t -> z(t, beta), x)
zpp(x, beta) = ForwardDiff.derivative(t -> zp(t, beta), x)

N = 100
beta = 0
xvec = LinRange(0, 1, N)
wvec = -zpp.(xvec, beta) .- 1

using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1], title="w", xlabel="x")
lines!(ax, xvec, wvec)
fig