using FiniteVolumes
using BenchmarkTools
using UnPack

include("../../src/numdiff/include_file.jl")

# Domain definition
Nx = 100
xmin, xmax = -4, 4
t0, tf = 0.0, 0.2
CFL_factor = 0.5
mesh = OneDMesh(Nx, xmin, xmax)
params = Parameters(mesh, t0, tf, CFL_factor)

function u0(x::Real)
    if x <= -1
        return zero(x)
    elseif x >=1
        return zero(x)
    else
        return one(x)
    end
end

u0(x::AbstractVector) = u0.(x)

aa = 3
equation = Equation(OneD(), 1, Scalar(), Advection(aa), u0)

# FV solution
# uinit = u0(mesh.x)
# dt = CFL_factor*mesh.dx/aa
# u = zero(uinit)
# for j in 1:Nx
#     u[j] = uinit[j] - aa*dt/mesh.dx*(uinit[j] - uinit[mod1(j-1,Nx)])
# end
# fluxvec = numflux(Euler(), Rusanov(), params, equation, cache, space_cache)
sol = solve(equation, params, Euler(), Rusanov(); maxiter=1, log_config=LogConfig(true, true, true, true, true));
# fnum = sol.log.fnumlog[end]
# plot(mesh.x, fnum, label="fnum")
# display(plot!(mesh.x, aa.*uinit, label="theo flux"))
uinit = sol.uinit
u = sol.u
dt = sol.dt

# using Plots
# plot(mesh.x, sol.log.fcontlog[end], label="fcont")
# scatter!(mesh.x, sol.u.*aa, label="u*a")
# #display(plot!(mesh.x, aa.*u, label="theo flux after"))
# display(plot!(mesh.x, aa.*uinit, label="theo flux"))


# Flux and entropy
Gexact = aa/2.0 * uinit.^2
etainit = uinit.^2/2.0
eta = u.^2/2.0

# Diffusion
D = zero(Gexact)

for j in 1:Nx
    D[j] = eta[j] - etainit[j] + dt/mesh.dx*(Gexact[j] - Gexact[mod1(j-1,Nx)])
end

plot(mesh.x, D, label="D")