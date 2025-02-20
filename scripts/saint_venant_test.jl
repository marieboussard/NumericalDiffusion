using FiniteVolumes
using Plots
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


sol = solve(equation, params, Euler(), Roe())#; log_config=LogConfig(true,true,true,true));

plot(mesh.x, znum, label="topo")
display(plot!(mesh.x, sol.u[:,1] .+ znum, label="water height"))

# plot(mesh.x, znum, label="topo")
# display(plot!(mesh.x, sol.uinit[:,1] .+ znum, label="water height"))

# display(plot(mesh.x, sol.u[:,2], label="water flow"))

# println("Burgers")
# @btime uinit = BurgersArticle.initcond(params.mesh.x)
# println("SaintVenant")
# @btime uinit = SaintVenantAtRest.initcond(params.mesh.x)



# println()
# println()
# println("Burgers")
# @show @allocated solve(BurgersArticle, params, Euler(), Rusanov(); maxiter=5);
# println()
# println("SaintVenant")
# @show @allocated solve(SaintVenantAtRest, params, Euler(), Rusanov(); maxiter=5);


# println()
# println()
# println("Burgers")
# @btime solve(BurgersArticle, params, Euler(), Rusanov(); maxiter=1000);
# println()
# println("SaintVenant")
# @btime solve(SaintVenantAtRest, params, Euler(), Rusanov(); maxiter=1000);

# println()
# println()
# println("Burgers")
# @allocated Integrator(BurgersArticle, params, Euler(), Rusanov(), 5, DefaultLogConfig);
# println()
# println("SaintVenant")
# @allocated Integrator(SaintVenantAtRest, params, Euler(), Rusanov(), 5, DefaultLogConfig);