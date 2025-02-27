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

# using Plots
# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="water height"))

# display(plot(mesh.x, sol.u[:,2], label="water flow"))

# function discretize_sourceterm2(::Pointwise, topo_source::TopoSource , v, mesh::Mesh, source_cache::TopoSourceCache)
#     @unpack Dznum = source_cache
#     @show @allocated s = similar(v)
#     @show @allocated nvar = ndims(Dznum)+1
#     # @show @allocated indices = [Colon() for _ in 1:nvar-1]
#     # FIRST EQUATION HAS ZERO SOURCE TERM
#     # @show @allocated s1 = view(s, indices..., 1)
#     @show @allocated s1 = selectdim(s, nvar, 1)
#     @show @allocated fill!(s1, zero(eltype(v)))
#     # for i in eachindex(s1)
#     #     s1[i] = zero(eltype(v))
#     # end
#     # OTHER EQUATIONS HAVE SPACE DERIVATED SOURCE TERMS
#     if nvar==2
#         for i in eachindex(Dznum)
#             s[i,2] = -v[i][1]*g*Dznum[i]
#         end
#     else
#         for r in 2:nvar
#             @show @allocated sr = selectdim(s, nvar, r)
#             for i in eachindex(s1)
#                 sr[i] = -g * v[i, 1] * Dznum[i][r-1]
#             end
#         end
#     end
#     s
# end

# integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig);

# @unpack equation, uinit, params, source_cache = integrator
# @unpack source = equation

# @show @allocated discretize_sourceterm2(source.source_discretize, source, uinit, params.mesh, source_cache)