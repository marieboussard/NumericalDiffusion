# f = -1

z0(x, beta) = x^2 / 2 + (beta - 1 / 2) * x + beta

Delta(beta) = beta^2 - 3 * beta + 1 / 4

x_plus(beta) = Delta(beta) >= 0 ? 1 / 2 - beta + sqrt(Delta(beta)) : 0
x_minus(beta) = Delta(beta) >= 0 ? 1 / 2 - beta - sqrt(Delta(beta)) : 0

beta_plus = (3 + 2sqrt(2)) / 2
beta_minus = (3 - 2sqrt(2)) / 2



betavec = LinRange(-5, 5, 10)

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1], title="Delta")
lines!(ax, betavec, Delta.(betavec))
lines!(ax, [beta_minus, beta_minus], [-1, 1])
lines!(ax, [beta_plus, beta_plus], [-1, 1])

fig

beta = 2.95
xvec = LinRange(0, 1, 100)

fig2 = Figure()
ax2 = Axis(fig2[1, 1], title="z0")
lines!(ax2, xvec, z0.(xvec, beta))
lines!(ax2, [x_plus(beta), x_plus(beta)], [-1, 1])
lines!(ax2, [x_minus(beta), x_minus(beta)], [-1, 1])

fig2

betavec = LinRange(-3, 3, 100)
#xvec = LinRange(0, 1, 100)

fig2 = Figure()
ax2 = Axis(fig2[1, 1], title="z0")
#lines!(ax2, xvec, z0.(xvec, beta))
lines!(ax2, betavec, x_minus.(betavec), label="x minus")
lines!(ax2, betavec, x_plus.(betavec), label="x plus")
axislegend()
fig2

