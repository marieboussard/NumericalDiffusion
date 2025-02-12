using FiniteVolumes
using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
xmin, xmax = 0, 1
t0, tf = 0.0, 10.0
CFL_factor = 0.5

mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

equation = SaintVenantAtRest
#equation = SaintVenantFlat

znum = z(equation.source, mesh.x)


# sol = solve(equation, params, Euler(), Rusanov(); maxiter=10);

# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="water height"))

# display(plot(mesh.x, sol.u[:,2], label="water flow"))

# println("Burgers")
# @btime uinit = BurgersArticle.initcond(params.mesh.x)
# println("SaintVenant")
# @btime uinit = SaintVenantAtRest.initcond(params.mesh.x)

# using FiniteVolumes: init_lake_at_rest
# using FiniteVolumes: zsinus

# function init_lake_at_rest2(x::AbstractVector, znum::AbstractVector; c=one(eltype(x)))
#     v = zeros(eltype(x), (length(x), 2))
#     for i in eachindex(x)
#         v[i,1] = max(zero(eltype(x)), c - znum[i])
#         v[i,2] = zero(eltype(x))
#     end
#     return v
# end


# zsinus2(x) = (-cos.(2*pi*freq * x) .+ 1)*height*0.5
# initcond2 = x -> init_lake_at_rest(x,zsinus);
# znum = zsinus.(mesh.x)
# initcond3 = x -> init_lake_at_rest2(x,znum)



# println()
# println()
# println("Burgers")
# @show @allocated solve(BurgersArticle, params, Euler(), Rusanov(); maxiter=5);
# println()
# println("SaintVenant")
# @show @allocated solve(SaintVenantAtRest, params, Euler(), Rusanov(); maxiter=5);


println()
println()
println("Burgers")
@btime solve(BurgersArticle, params, Euler(), Rusanov(); maxiter=1000);
println()
println("SaintVenant")
@btime solve(SaintVenantAtRest, params, Euler(), Rusanov(); maxiter=1000);

# println()
# println()
# println("Burgers")
# @allocated Integrator(BurgersArticle, params, Euler(), Rusanov(), 5, DefaultLogConfig);
# println()
# println("SaintVenant")
# @allocated Integrator(SaintVenantAtRest, params, Euler(), Rusanov(), 5, DefaultLogConfig);


# @show @allocated integrator = Integrator(equation, params, Euler(), Rusanov(), 100, DefaultLogConfig);

# using FiniteVolumes: view_stencil!

# view_stencil!(integrator, 2)
# #VSCodeServer.@profview solve(equation, params, Euler(), Rusanov())
# const g=9.8
# function discretize_sourceterm2!(::Pointwise, integrator::Integrator)
#     @unpack uprev, cache, source_cache = integrator
#     @views h = uprev[:,1]
#     for i in eachindex(h)
#         cache.sourceterm[i,2] = -h[i]*g*source_cache.Dznum[i]
#     end
# end

# function numflux2!(integrator::Integrator)
#     for i ∈ 1:integrator.params.mesh.Nx
#         view_stencil!(integrator, i)
#         numflux!(integrator.time_scheme, integrator, i+1)
#     end
#     for j in 1:integrator.equation.p
#         integrator.fnum[end,j] = integrator.fnum[1,j]
#     end
#     # @show @allocated integrator.fnum[end,:] .= integrator.fnum[1,:]
#     nothing
# end


# function loopfooter2!(integrator::Integrator)
#     @unpack cache, equation, uprev, u, fcont = integrator
#     @unpack source = equation

#     @show @allocated integrator.t += integrator.dt
#     @show @allocated integrator.niter += 1
#     @show @allocated uprev .= u
#     @show @allocated fcont .= flux(equation.funcs, u)
#     # STORING FLUX DERIVATIVE IN SCALAR CASE
#     if equation.p ==1
#         @show @allocated cache.Dfcont .= Dflux(equation.funcs, u)
#     end
#     # UPDATING SOURCE TERM
#     if has_source(source)
#         @show @allocated discretize_sourceterm!(source.source_discretize, integrator)
#     end
#     @show @allocated update_log!(integrator)
# end