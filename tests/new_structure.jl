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

function numflux5!(::Rusanov, integrator::Integrator, i, args...)
    @show @allocated view_stencil!(integrator, i-1)
    @unpack equation, cache, fnum, fcont, uprev = integrator
    @unpack stencil = integrator.cache
    # @show @allocated uL   = view(uprev, stencil[1], :)
    # @show @allocated uR   = view(uprev, stencil[2], :)

    # @show @allocated @views fL = fcont[stencil[1], :]
    # @show @allocated @inbounds fL   = view(fcont, stencil[1], :)
    # @show typeof(fL)
    # @show @allocated fR   = view(fcont, stencil[2], :)
    
    # @show @allocated fnum_i = view(fnum, i, :)
    # @show @allocated fnum_i .= (fL .+ fR) ./ 2 - cache.cfl_loc./ 2 * (uR .- uL)
    @show @allocated fnum[i] = (fcont[stencil[1]] + fcont[stencil[2]]) *0.5

    for j in 1:equation.p
        @show @allocated fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5
        @show @allocated fnum[i,j] -= cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
    end

    @inbounds for j in 1:equation.p
        @show @allocated uL   = view(uprev, stencil[1], j)
        @show @allocated uR   = view(uprev, stencil[2], j)

        @show @allocated @views fL = fcont[stencil[1], j]
        @show @allocated @inbounds fL   = view(fcont, stencil[1], j)
        @show typeof(fL)
        @show @allocated fR   = view(fcont, stencil[2], j)
        

        @show @allocated fnum_ij = view(fnum, i, j)
        @show typeof(fnum_ij)
        @show @allocated @. fnum_ij = (fL + fR) *0.5
        @show @allocated @. fnum_ij -= cache.cfl_loc*0.5 * (uR - uL)
    end

end