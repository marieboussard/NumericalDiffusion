using Plots

include("../src/tools/equation.jl")
include("../src/tools/method.jl")
include("../src/tools/domain.jl")
include("../src/tools/fv_solution.jl")

Nx = 10

omega = createUnitInterval(Nx, 0.0, 0.4)
#v0 = v0_lake_at_rest(omega.x, Bump_zb())
# plot(omega.x, zb(Bump_zb(), omega.x), label="zb")
# plot!(omega.x, [e[1] for e in v0] .+ zb(Bump_zb(), omega.x), label="h+zb")
# xlabel!("x")

CFL_number = 0.5
sol = fv_solve(omega, x -> v0_lake_at_rest(x, Bump_zb()), SaintVenant(), Rusanov(CFL_number))

display(plot(omega.x, [sol.u_approx[end-1][i][1] for i in 1:Nx], label="h"))
plot(omega.x, [sol.u_approx[end-1][i][2] for i in 1:Nx], label="hu")