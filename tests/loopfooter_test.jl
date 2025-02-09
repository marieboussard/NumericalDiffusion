using FiniteVolumes
# using Plots
using BenchmarkTools
using UnPack
using Profile

Profile.init(delay = 0.0001)

# Domain definition
Nx = 1000
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Burgers equation
f(u) = u.^2/2
Df(u) = u
u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
equation = Equation(1, f, Df, u0)

# 16 allocations
integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig)

# 2 allocations
function loopfooter2!(integrator::Integrator)
    integrator.t += integrator.dt
    integrator.niter += 1
    integrator.uprev .= integrator.u
    integrator.fcont .= integrator.equation.funcs.flux.(integrator.u)
    update_log!(integrator)
end

# 6 allocations
function loopheader2!(integrator::Integrator)
    dt_CFL!(integrator)
end

function performstep2!(integrator::Integrator)
    @unpack dx = integrator.params.mesh
    @unpack u, uprev, dt, fnum = integrator
    numflux!(integrator)
    @views fluxforward = fnum[2:end,:]
    @views fluxbackward = fnum[1:end-1,:]
    @. u = uprev - dt / dx * (fluxforward .- fluxbackward)
end