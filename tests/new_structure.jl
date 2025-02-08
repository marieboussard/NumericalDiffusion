using FiniteVolumes
# using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 10000
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
sol = solve(equation, params, Euler(), Rusanov())

# plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
# plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.params.tf))


integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig)

#u, fu = FiniteVolumes.view_stencil!(integrator, 3)

function numflux2!(::Rusanov, integrator::Integrator, u, i, args...)
    @unpack equation, cache, fnum = integrator
    @unpack flux = equation

    uL = view(u, 1, :)
    uR = view(u, 2, :)

    cache.cfl_loc = maximum(abs, u)
    
    fnum_i = view(fnum, i, :)

    @show @allocated @. fnum_i = (flux(uL) + flux(uR)) / 2 
    
    @show @allocated @. fnum_i -= cache.cfl_loc/ 2 * (uR - uL)
end

function numflux3!(::Rusanov, integrator::Integrator, u, i, args...)
    @unpack equation, cache, fnum = integrator
    @unpack flux = equation

    uL = view(u, 1, :)
    uR = view(u, 2, :)

    cache.cfl_loc = maximum(abs, u)
    
    @show @allocated fnum[i,:] .= (flux(uL) .+ flux(uR)) ./ 2 
    
    @show @allocated @. fnum[i,:] .-= cache.cfl_loc./ 2 * (uR .- uL)
end

using FiniteVolumes: view_stencil!
function numflux4!(::Rusanov, integrator::Integrator, i)
    @unpack equation, cache, fnum, fcont, uprev = integrator
    @unpack stencil = integrator.cache
    @show @allocated view_stencil!(integrator, i)
    @show @allocated cache.cfl_loc = maximum([abs.(uprev[stencil[1]]), abs.(uprev[stencil[2]])])

    @show @allocated @inbounds fnum[i] = uprev[stencil[1]]
    @show @allocated @inbounds fnum[i] -= uprev[stencil[2]]
    @show @allocated @inbounds fnum[i] +=  (fcont[stencil[1]] + fcont[stencil[2]]) 
    @show @allocated fnum[i] *= 0.5
    nothing
end

