using NumericalDiffusion

u0(x::Real) = x<=0 ? 1.0 : -1.0
u0(x::AbstractVector) = u0.(x)

Nx = 100
xmin, xmax = -2.0, 2.0
t0, tf = 0.0, 0.4 
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

equation = Equation(OneD(), 1, System(), SaintVenant(), u0)
time_scheme = Euler()
space_scheme = Rusanov()

sol = solve(equation, params, time_scheme, space_scheme)