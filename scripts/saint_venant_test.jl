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

equation = SaintVenantAtRest
#equation = SaintVenantFlat

znum = z(equation.source, mesh.x)


sol = solve(equation, params, Euler(), Roe());#; log_config=LogConfig(true,true,true,true));

using Plots
plot(mesh.x, znum, label="topo")
display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="water height"))

# display(plot(mesh.x, sol.u[:,2], label="water flow"))
# const nvar = 2

# function discretize_sourceterm2!(::Pointwise, integrator::Integrator)
#     @unpack u, cache, source_cache = integrator
#     @unpack Dznum, nvar = source_cache
#     @unpack sourceterm = cache
#     # @show @allocated nvar = ndims(znum)+1
#     @show @allocated indices = ntuple(i -> Colon(), nvar - 1)
#     #@show @allocated indices = indices = [Colon() for i in 1:nvar-1]
#     # FIRST EQUATION HAS ZERO SOURCE TERM
#     @show @allocated h = view(u, indices..., 1)
#     @show @allocated h = selectdim(u, nvar, 1)
#     @show @allocated s1 = view(sourceterm, indices..., 1)
#     @show @allocated s1 = selectdim(sourceterm, nvar, 1)
#     @show @allocated fill!(s1, zero(eltype(u)))
#     # OTHER EQUATIONS HAVE SPACE DERIVATED SOURCE TERMS
#     if nvar==2
#         for i in eachindex(Dznum)
#             @show @allocated sourceterm[i,2] = -g*h[i]*Dznum[i]
#         end
#     else
#         for r in 2:nvar
#             @show @allocated sr = selectdim(sourceterm, nvar, r)
#             for i in eachindex(s1)
#                 @show @allocated sr[i] = -g * h[i] * Dznum[i][r-1]
#             end
#         end
#     end
# end
# integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig);

# @show @allocated discretize_sourceterm2!(integrator.equation.source.source_discretize, integrator)

# @unpack equation, uinit, params, source_cache = integrator
# @unpack source = equation

# @show @allocated discretize_sourceterm2(source.source_discretize, source, uinit, params.mesh, source_cache)