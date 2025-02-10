using FiniteVolumes
using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
xmin, xmax = -2, 2
t0, tf = 0.0, 0.4
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

# Burgers equation
f(u) = u.^2/2
Df(u) = u
u0(x::Real) = x <= 0 ? -2 -x : 3 - 3 / 2 * x
u0(x::AbstractVector) = u0.(x)
equation = Equation(1, Scalar(), EquationFun(f, Df), u0)
sol = solve(equation, params, Euler(), Rusanov(); log_config=LogConfig(true, true, true));

plt = plot(sol.params.mesh.x, sol.uinit, label=string(sol.params.t0))
display(plot!(plt, sol.params.mesh.x, sol.u, label=string(sol.params.tf)))


integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig)

# #u, fu = FiniteVolumes.view_stencil!(integrator, 3)


using FiniteVolumes: view_stencil!


# function numflux2!(::Rusanov, integrator::Integrator, i, args...)
#     @unpack equation, cache, fnum, fcont, Dfcont, uprev = integrator
#     @unpack stencil = integrator.cache
#     view_stencil!(integrator, i-1)
#     for j in 1:equation.p
#         cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))
#         fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
#     end
# end

# function numflux3!(::Rusanov, integrator::Integrator, i, args...)
#     @unpack equation, cache, fnum, fcont, Dfcont, uprev = integrator
#     @unpack stencil = integrator.cache
#     view_stencil!(integrator, i-1)
#     for j in 1:equation.p
#         CFL_local!(equation.eqtype, integrator, j)
#         fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5 - cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
#     end
# end








# # function numflux4!(::Rusanov, integrator::Integrator, i)
# #     @unpack equation, cache, fnum, fcont, uprev = integrator
# #     @unpack stencil = integrator.cache
# #     @show @allocated view_stencil!(integrator, i)
# #     @show @allocated cache.cfl_loc = maximum([abs.(uprev[stencil[1]]), abs.(uprev[stencil[2]])])

# #     @show @allocated @inbounds fnum[i] = uprev[stencil[1]]
# #     @show @allocated @inbounds fnum[i] -= uprev[stencil[2]]
# #     @show @allocated @inbounds fnum[i] +=  (fcont[stencil[1]] + fcont[stencil[2]]) 
# #     @show @allocated fnum[i] *= 0.5
# #     nothing
# # end

# function numflux5!(::Rusanov, integrator::Integrator, i, args...)
#     @show @allocated view_stencil!(integrator, i-1)
#     @show @allocated @unpack equation, cache, fnum, fcont, Dfcont, uprev = integrator
#     @show @allocated @unpack stencil = integrator.cache
#     # @show @allocated uL   = view(uprev, stencil[1], :)
#     # @show @allocated uR   = view(uprev, stencil[2], :)

#     # @show @allocated @views fL = fcont[stencil[1], :]
#     # @show @allocated @inbounds fL   = view(fcont, stencil[1], :)
#     # @show typeof(fL)
#     # @show @allocated fR   = view(fcont, stencil[2], :)
    
#     # @show @allocated fnum_i = view(fnum, i, :)
#     # @show @allocated fnum_i .= (fL .+ fR) ./ 2 - cache.cfl_loc./ 2 * (uR .- uL)
#     #@show @allocated fnum[i] = (fcont[stencil[1]] + fcont[stencil[2]]) *0.5

#     for j in 1:equation.p
#         @show @allocated cache.cfl_loc = max(abs(Dfcont[stencil[1],j]), abs(Dfcont[stencil[2],j]))

#         @show @allocated fnum[i,j] = (fcont[stencil[1], j] + fcont[stencil[2],j]) *0.5
#         @show @allocated fnum[i,j] -= cache.cfl_loc.*0.5 * (uprev[stencil[2],j] - uprev[stencil[1],j])
#     end

#     # @inbounds for j in 1:equation.p



#     #     @show @allocated uL   = view(uprev, stencil[1], j)
#     #     @show @allocated uR   = view(uprev, stencil[2], j)

#     #     @show @allocated @views fL = fcont[stencil[1], j]
#     #     @show @allocated @inbounds fL   = view(fcont, stencil[1], j)
#     #     @show typeof(fL)
#     #     @show @allocated fR   = view(fcont, stencil[2], j)
        

#     #     @show @allocated fnum_ij = view(fnum, i, j)
#     #     @show typeof(fnum_ij)
#     #     @show @allocated @. fnum_ij = (fL + fR) *0.5
#     #     @show @allocated @. fnum_ij -= cache.cfl_loc*0.5 * (uR - uL)
#     # end

# end

# # function performstep2!(integrator::Integrator)
# #     @unpack dx = integrator.params.mesh
# #     @unpack u, uprev, dt, fnum, equation, params = integrator
# #     @show @allocated numflux!(integrator)
# #     for i in 1:params.mesh.Nx
# #         for j in 1:equation.p
# #             @show @allocated u[i,j] = uprev[i,j] - dt / dx * (fnum[i+1,j] - fnum[i,j])
# #         end
# #     end
# # end