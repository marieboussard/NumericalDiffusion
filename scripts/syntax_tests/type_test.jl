using FiniteVolumes
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
xmin, xmax = 0, 1
t0, tf = 0.0, 0.5
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# SINUSOIDAL TOPOGRAPHY
const freq = 1.0
const height = 0.5
z(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
Dz(x) = pi*freq*(sin(2*pi*freq * x))*height

equation = saintvenant_with_topo(z, Dz; sourcedisc=HRDisc())
scheme = HR(Rusanov())

integrator = Integrator(equation, params, Euler(), scheme, 100, DefaultLogConfig)

function numtest!(scheme::HR, integrator::Integrator, j::Int, args...)
    @unpack equation, params, space_cache, uprev = integrator
    @unpack Nx = params.mesh
    hplus!(integrator, j)
    hminus!(integrator, j)
    @unpack hplus, hminus, uplus, uminus, fplus, fminus = space_cache
    uplus[1], uplus[2] = hplus[j], hplus[j]*uprev[mod1(j+1,Nx),2]
    uminus[1], uminus[2] = hminus[j], hminus[j]*uprev[j,2]
    nothing
end

function numtest2!(scheme::HR, integrator::Integrator, j::Int, args...)
    @unpack equation, params, space_cache, uprev = integrator
    @unpack Nx = params.mesh
    hplus!(integrator, j)
    hminus!(integrator, j)
    @unpack hplus, hminus, uplus, uminus, fplus, fminus = space_cache
    uplus .= hplus[j] .*(one(eltype(hplus)), uprev[mod1(j+1,Nx),2])
    uminus .= hminus[j] .*(one(eltype(hplus)), uprev[j,2])
    nothing
end