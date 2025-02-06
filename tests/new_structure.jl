# Domain definition

Nx = 10
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Burgers equation
f(u) = u.^2/2
u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x

sol = solve(Equation(1, f, u0), params, Euler(), Rusanov())