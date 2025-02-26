using FiniteVolumes
# using Plots
using BenchmarkTools
using UnPack

# Domain definition
Nx = 100
Ny = 100
xmin, xmax, ymin, ymax = -2, 2,-2,2
t0, tf = 0.0, 2.0
CFL_factor = 0.5

mesh = TwoDCartesian(Nx, Ny, xmin, xmax, ymin, ymax)
params = Parameters(mesh, t0, tf, CFL_factor)

using FiniteVolumes: SaintVenant2Flat
equation = SaintVenant2Flat

# sol = solve(equation, params, Euler(), Rusanov2D(); maxiter=1);

# display(heatmap(mesh.x, mesh.y, sol.u[:,:,1]))
# display(heatmap(mesh.x, mesh.y, sol.u[:,:,2]))
# heatmap(mesh.x, mesh.y, sol.u[:,:,3])

integrator = Integrator(equation, params, Euler(), Rusanov2D(), 1, DefaultLogConfig);

# function loopfooter2!(integrator::Integrator)
#     @unpack equation, uprev, u = integrator
#     @unpack source = equation
#     @show @allocated integrator.t += integrator.dt
#     @show @allocated integrator.niter += 1
#     @show @allocated uprev .= u
#     @show @allocated update_flux!(equation.dim, integrator)
#     @show @allocated update_cflcache!(equation.dim, equation.eqtype, equation.funcs, integrator)
#     # UPDATING SOURCE TERM
#     if has_source(source)
#         discretize_sourceterm!(source.source_discretize, integrator)
#     end
#     @show @allocated update_log!(integrator)
#     nothing
# end

# function update_flux2!(::TwoD, integrator::Integrator)
#     @unpack equation, u = integrator
#     @unpack fcont, hcont = integrator.fcont
#     @show @allocated fcont .= flux_f(equation.funcs, u)
#     @show @allocated hcont .= flux_h(equation.funcs, u)
# end

# function flux_f2(::SaintVenant2D, v)
#     res = similar(v)
#     g_half = g * 0.5 
#     h = view(v, :, :, 1)
#     hu = view(v, :, :, 2)
#     hv = view(v, :, :, 3)
#     #for i in eachindex(h)
#     for I in CartesianIndices(h)
#         if h[I] > treshold
#             # hu = v[i, 2]
#             # hv = v[i, 3]
#             res[I, 1] = hu[I]
#             res[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
#             res[I, 3] = hu[I]*hv[I] / h[I]
#         else
#             res[I, 1] = zero(eltype(v))
#             res[I, 2] = zero(eltype(v))
#             res[I, 3] = zero(eltype(v))
#         end
#     end
#     res
# end

# function flux_f2!(::SaintVenant2D, integrator)
#     @unpack fcont, hcont = integrator.fcont
#     @unpack u = integrator
#     g_half = g * 0.5 
#     h = view(u, :, :, 1)
#     hu = view(u, :, :, 2)
#     hv = view(u, :, :, 3)
#     for I in CartesianIndices(h)
#         if h[I] > treshold
#             fcont[I, 1] = hu[I]
#             fcont[I, 2] = hu[I]^2 / h[I] + g_half * h[I]^2
#             fcont[I, 3] = hu[I]*hv[I] / h[I]
#         else
#             fcont[I, 1] = zero(eltype(u))
#             fcont[I, 2] = zero(eltype(u))
#             fcont[I, 3] = zero(eltype(u))
#         end
#     end
# end